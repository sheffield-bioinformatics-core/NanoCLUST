#!/usr/bin/env python

"""Create results report."""

import argparse
from concurrent.futures import process
from curses import meta
import uuid
from aplanat import report
import pandas as pd
import os

#need our own UoSReport class

class UoSReport(report.HTMLReport):
    """Report template for Sheffield Bioinformatics NanoCLUST workflow."""

    def __init__(
            self, title, workflow, commit=None, revision=None,
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
            title=title, lead=lead, require_keys=require_keys, style=style)
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
    info=pd.read_excel(file, usecols=range(0,5))
    #return row with a barcode as a Series
    relevant_row=info.loc[info['Barcode'] == barcode].transpose()
    #move row names into a column
    print(barcode)
    relevant_row.index.name = 'Metadata'
    relevant_row.reset_index(inplace=True)
    #rename columns
    relevant_row.columns=['Metadata', 'Sample Information']
    return relevant_row

def read_abundance_results(file):
    abundance_results=pd.read_csv(file)
    abundance_results.columns=['Detected Species', 'Relative Abundance (%)', 'Number of Reads']
    return abundance_results

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
        "--output", required=True,
        help="Report output file name")
    parser.add_argument(
        "--barcode", required=True,
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
    args = parser.parse_args()

    metadata_table=read_patient_info(args.info, args.barcode)
    title="Patient " + metadata_table.iloc[0,1] + " Report"

    results_table=read_abundance_results(args.infile)
    total_reads=results_table['Number of Reads'].sum()

    positive,negative=process_controls(args.controls)

    report = UoSReport(
        title, "NanoCLUST",
        revision=args.revision, commit=args.commit)

    section=report.add_section()

    section.markdown('''
    ### Results

    Total reads in this sample: {0}
    
    Matched species:
    '''.format(total_reads))

    section.table(results_table)
    
    if metadata_table.loc[metadata_table['Metadata'] == 'Infection type', 'Sample Information'].iloc[0] == 'bacterial':
        database_info="16s bacterial sequencing results were compared against 16S & 18S database, build 18 Jan 2022."
    else:
        database_info="ITS2 fungal sequencing results were compared against ITS2 database, build 15 Mar 2022."

    section.markdown('''
    Sample was sequenced on a ONT GridION Mk1. 
    Sequencing data was processed and analysed using a custom nanoclust pipeline. 
    {0}
    '''.format(database_info))

    section=report.add_section()
    section.markdown('''
    ### Run QC


    **NEGATIVE CONTROL**
    ''')

    if negative is not None:
        section.markdown('''
        Total reads in negative control: {0} 
        '''.format(negative['Number of Reads'].sum()))

        section.table(negative)
    else:
        section.markdown('''
        Negative control was not provided.
        ''')

    section.markdown('''
    **POSITIVE CONTROL**
    ''')

    if positive is not None:
        section.markdown('''
        Total reads in positive control: {0}
        '''.format(positive['Number of Reads'].sum()))
        
        section.table(positive)
    else:
        section.markdown('''
        Positive control was not provided.
        ''')

    section.markdown('''
    comment: green/amber/reds 
    ''')

    section=report.add_section()
    section.markdown('''
    ### Metadata

    This section displays the basic metadata.
    ''')

    section.table(metadata_table)

    section=report.add_section()
    run_id='example run ID'
    barcoding_kit='example barcoding kit'
    demux_method=args.demux
    species_database=metadata_table['Sample Information'].iloc[4]
    clustering_size=args.clustering_size
    section.markdown('''
    ### Run parameters

    **Run ID**: (should be able to pull this out from the run report) {0}

    **Barcoding kit**: (should be able to pull this out from the run report; report_RUNID.pdf/Run Parameters/Kit) {1}

    **Demultiplex method**: {2}

    **Species Database**: {3}

    **Clustering Size**: {4}
    '''.format(run_id, barcoding_kit, demux_method, species_database, clustering_size))

    #write report
    report.write(args.output)


if __name__ == "__main__":
    main()