#!/usr/bin/env python

"""Create results report."""

import argparse
from concurrent.futures import process
from curses import meta
import uuid
from aplanat import report
import pandas as pd
import os
import numpy as np
from jinja2 import Template
from collections import OrderedDict
import pkg_resources
from bokeh.resources import INLINE

class HTMLReport(report.HTMLSection):
    """Generate HTML Report from a series of bokeh figures.
    Items added to the report take an optional key argument, adding items
    with the same key allows an update in place whilst maintaining the order
    in which items were added. Items can be grouped into sections for easier
    out of order addition.
    """

    def __init__(self, title="", lead="", report_template="", require_keys=False, style='ont'):
        """Initialize the report item collection.
        :param title: report title.
        :param lead: report strapline, shown below title.
        :param require_keys: require keys when adding items.
        """
        super().__init__(require_keys=require_keys)
        self.title = title
        self.lead = lead
        self.sections = OrderedDict()
        self.sections['main'] = self
        self.style = style

        self.CSS_RESOURCES = dict(
                bootstrap='bootstrap.min.css',
                datatables='simple-datatables_latest.css',
                ont='custom-ont.css',
                ond='custom-ond.css',
                epi2me='custom-epi2me.css'
        )

        template = report_template
        with open(template, 'r', encoding="UTF-8") as fh:
            template = fh.read()
        self.template = Template(template)

    def add_section(self, key=None, section=None, require_keys=False):
        """Add a section (grouping of items) to the report.
        :param key: unique key for section.
        :param section: `HTMLSection` to add rather than creating anew.
        :returns: the report section.
        """
        if key is None:
            key = str(uuid.uuid4())
        section = report._maybe_new_report(section, require_keys=require_keys)
        self.sections[key] = section
        return self.sections[key]

    def render(self):
        """Generate HTML report containing figures."""
        bokeh_resources = INLINE.render()

        libs = []
        css_files = ['bootstrap', 'datatables', self.style]
        CSS_RESOURCES = [self.CSS_RESOURCES[file] for file in css_files]
        JS_RESOURCES = ['simple-datatables_latest.js']
        for resources, stub in (
                [CSS_RESOURCES, "<style>{}</style>"],
                [JS_RESOURCES, '<script type="text/javascript">{}</script>']):
            for res in resources:
                fn = pkg_resources.resource_filename(
                    "aplanat.report", 'data/{}'.format(res))
                with open(fn, encoding="UTF-8") as fh:
                    libs.append(stub.format(fh.read()))

        all_scripts = list()
        all_divs = list()
        for sec_name, section in self.sections.items():
            scripts, divs = section.components()
            all_scripts.extend(scripts)
            all_divs.extend(divs)

        script = '\n'.join(all_scripts)
        divs = '\n'.join(all_divs)
        return self.template.render(
            title=self.title, lead=self.lead, logo="UoS",
            bokeh_resources=bokeh_resources, resources="\n".join(libs),
            script=script, div=divs)

    def write(self, path):
        """Write html report to file."""
        with open(path, "w", encoding='utf8') as outfile:
            outfile.write(self.render())

#need our own UoSReport class

class UoSReport(HTMLReport):
    """Report template for Sheffield Bioinformatics NanoCLUST workflow."""

    def __init__(
            self, title, report_template, workflow=None, commit=None, revision=None,
            require_keys=False, about=True, style='ont'):
        """Initialize the report item collection.
        :param workflow: workflow name (NanoCLUST)
        :param title: report title.
        :param require_keys: require keys when adding items.
        """
        if commit is None:
            commit = "unknown"
        if revision is None:
            revision = "unknown"
        self.commit = commit
        self.revision = revision
        self.workflow = workflow
        self.about = about
        self.style = style

        lead = (
            "Results generated through the {} Nextflow workflow "
            "provided by University of Sheffield Bioinformatics.".format(workflow))
        super().__init__(
            title=title, lead=lead, report_template=report_template, require_keys=require_keys, style=style)
        self.tail_key = str(uuid.uuid4())

    def render(self):
        """Generate HTML report containing figures."""
        # delete and re-add the tail (in case we are called twice,
        # and something was added).
        try:
            del self.sections[self.tail_key]
        except KeyError:
            pass

        if self.about:
            self.add_section(key=self.tail_key).markdown("""
    ### About
    This report was produced using the
    [sheffield-bioinformatics-core/{0}](https://github.com/sheffield-bioinformatics-core/{0}).  The
    workflow can be run using `nextflow sheffield-bioinformatics-core/{0} --help`
    **Version details** *Revision*: {1} *Git Commit*: {2}
    ---
    """.format(self.workflow, self.revision, self.commit))
        return super().render()

def read_patient_info(file, barcode):
    #funtion parsing CSV patient file and looking up info for relevant barcode
    info=pd.read_excel(file, usecols=range(0,7))
    #return row with a barcode as a Series
    if barcode=="discontinued":
        relevant_row=[]
        relevant_rows=info.loc[info['Status'] == barcode]
        for index,row in relevant_rows.iterrows():
            single_row=row.transpose()
            single_df=single_row.reset_index()
            single_df.columns=['Metadata', 'Sample Information']
            relevant_row.append(single_df)
    else:
        relevant_rows=info.loc[info['Barcode'] == barcode]
        #move row names into a column
        relevant_row=relevant_rows.transpose()
        relevant_row.index.name = 'Metadata'
        relevant_row.reset_index(inplace=True)
        #rename columns
        relevant_row.columns=['Metadata', 'Sample Information']
    return relevant_row

def read_abundance_results(file):
    abundance_results=pd.read_csv(file)
    abundance_results.columns=['Detected Species', 'Relative Abundance (%)', 'Number of Reads']
    abundance_results['Relative Abundance (%)']=abundance_results['Relative Abundance (%)'].apply(np.around)
    abundance_results_t3=abundance_results.head(n=3)

    return abundance_results_t3

def process_controls(controls):
    for i in controls:
        if "positive" in i:
            if os.stat(i).st_size == 0:
                positive = None
            else:
                positive = read_abundance_results(i)
        else:
            if os.stat(i).st_size == 0:
                negative = None
            else:
                negative = read_abundance_results(i)

    return positive,negative

def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile", default='unknown',
        help='Table file with classification and abundance results')
    parser.add_argument(
        "--output", default='unknown',
        help="Report output file name")
    parser.add_argument(
        "--barcode", default='discontinued',
        help="barcode identifier for the patient sample")
    parser.add_argument(
        "--info", required=True,
        help="Experiment information file mapping patient ID with a barcode")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--demux", default='unknown',
        help='demultiplexing method and software version')
    parser.add_argument(
        "--clustering_size", default='unknown',
        help="Amount of reads used for UMAP HDBSCAN clustering")
    parser.add_argument(
        "--controls", default='unknown', nargs=2,
        help="File names for positive and negative control results")
    parser.add_argument(
        "--reads_count", default='0',
        help="Reads count after quality control")
    parser.add_argument(
        "--kit", default='unknown',
        help="Kit used for barcoding and demultiplexing")
    parser.add_argument(
        "--report_template",
        help="path to report template")
    args = parser.parse_args()
    if args.barcode=="discontinued":

        metadata_table_list=read_patient_info(args.info, args.barcode)
        for patient in metadata_table_list:
            title="Patient " + patient.iloc[0,1] + " Report"
            report = UoSReport(
                title=title, report_template=args.report_template, about=False)

            section=report.add_section()
            section.markdown('''
            ### Sample Information

            This section displays the basic metadata.
            ''')

            section.table(patient.iloc[[0,1,2,4,6]])

            section=report.add_section()

            assay_type=patient.iloc[4]['Sample Information']
            if assay_type == '16S':
                assay_info = 'Bacterial 16s'
            else:
                assay_info = 'Fungal ITS2'

            section.markdown('''
            ### Results

            <font color="red">**{0} rRNA NOT detected**</font>
            '''.format(assay_info))

            report.write("patient_report_" + str(patient.iloc[1,1]) + ".html")
    else:        
        metadata_table=read_patient_info(args.info, args.barcode)
        title="Patient " + metadata_table.iloc[0,1] + " Report"

        if args.infile == "input.1":
            results_table = None
        else:
            results_table=read_abundance_results(args.infile)

        positive,negative=process_controls(args.controls)

        report = UoSReport(
            title=title, workflow="NanoCLUST", report_template=args.report_template,
            revision=args.revision, commit=args.commit)

        section=report.add_section()
        section.markdown('''
        ### Sample Information

        This section displays the basic metadata.
        ''')

        section.table(metadata_table.iloc[[0,2,4,6]])

        section=report.add_section()

        section.markdown('''
        ### Results

        Total reads in this sample: {0}
        
        Matched species:
        '''.format(args.reads_count))

        assay_type=metadata_table.loc[metadata_table['Metadata'] == 'Assay', 'Sample Information'].iloc[0]
        if assay_type == '16S':
            database_info="16s bacterial sequencing results were compared against 16S & 18S database, build 18 Jan 2022."
            infection_type='bacterial'
            assay_info='Bacterial 16s'
        else:
            database_info="ITS2 fungal sequencing results were compared against ITS2 database, build 15 Mar 2022."
            infection_type='fungal'
            assay_info='Fungal ITS2'

        if results_table is not None:
            section.table(results_table)
        else:
            section.markdown('''
            <font color="red">**{0} rRNA NOT detected**</font>
            '''.format(assay_info))
        

        section=report.add_section()
        section.markdown('''
        ### Run QC


        **NEGATIVE CONTROL**
        ''')
        comment=0
        if negative is not None:
            comment+=10
            section.markdown('''
            Total reads in negative control: {0} 
            '''.format(negative['Number of Reads'].sum()))

            section.table(negative)
        else:
            section.markdown('''
            No species detected in negative control.
            ''')

        section.markdown('''
        **POSITIVE CONTROL**
        ''')

        if positive is not None:
            comment+=1
            section.markdown('''
            Total reads in positive control: {0}
            '''.format(positive['Number of Reads'].sum()))
            
            section.table(positive)
        else:
            comment+=2
            section.markdown('''
            No species detected in positive control.
            ''')

        if comment == 1:
            section.markdown('''
            comment: <font color="green">QC for this sample was **successful**</font>
            ''')
        elif comment == 11:
            section.markdown('''
            comment: <font color="orange">QC for this sample shows the presence of {0} reads in the negative control. Check if there is overlap with any detected pathogen in the sample that may indicate contamination.</font>
            '''.format(infection_type))
        else:
            section.markdown('''
            comment: <font color="red">QC for this sample has **failed** and results cannot be validated</font>
            ''')

        section=report.add_section()
        run_id='example run ID'
        barcoding_kit=args.kit
        print(barcoding_kit)
        demux_method=args.demux
        species_database=metadata_table['Sample Information'].iloc[4]
        clustering_size=args.clustering_size
        section.markdown('''
        ### Run parameters

        **Run ID**: (should be able to pull this out from the run report) {0}

        **Barcoding kit**: {1}

        **Demultiplex method**: {2}

        **Species Database**: {3}

        **Clustering Size**: {4}

        **Sample barcode**: {5}

        Sample was sequenced on a ONT GridION Mk1. 
        Sequencing data was processed and analysed using a custom nanoclust pipeline.
        {6}

        '''.format(run_id, barcoding_kit, demux_method, species_database, clustering_size, metadata_table.iloc[1,1], database_info))

        #write report
        report.write(args.output + "_" + str(metadata_table.iloc[1,1]) + ".html")


if __name__ == "__main__":
    main()