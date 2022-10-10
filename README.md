# NanoCLUST

**De novo clustering and consensus building for ONT 16S sequencing data**.

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install [`conda`](https://conda.io/miniconda.html) or [`docker`](https://docs.docker.com/engine/installation/)

iii. Clone the NanoCLUST repository and test the pipeline on a minimal dataset with a single command and conda/docker profiles (docker profiles currently outdated, so using **conda** is recommended)

### Download appropriate seqmatch or kraken2 databases needed for your analysis

To build a bacterial 16S kraken2 database run:
```bash
wget ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
gzip -d bacteria.16SrRNA.fna.gz
tar -xf new_taxdump.tar.gz
kraken2-build --download-taxonomy --db 16S_NCBI_30_08_2022
kraken2-build --add-to-library bacteria.16SrRNA.fna --db 16S_NCBI_30_08_2022
kraken2-build --build --db 16S_NCBI_30_08_2022
```
To run BLAST 16S rRNA classification download appropriate databases:

```bash
mkdir db db/taxdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && tar -xzvf 16S_ribosomal_RNA.tar.gz -C db
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzvf taxdb.tar.gz -C db/taxdb
```
### Prepare your custom config profile

```
profiles {
    conda {
        conda.cacheDir = "/path/to/conda/envs"
    }
}

params {
    classification = "kraken2"
    db = "/path/to/kraken2/database"
    tax = "/path/to/new_taxdump/rankedlineage.dmp"
    email = "your.email@email.com"
}
```

### Run a test

```bash
#Using conda profile (recommended).
nextflow run main.nf -profile test,conda
```

iv. Start running your own analysis!

Run a single sample analysis (default requires 32G RAM):

```bash
nextflow run /path/to/NanoCLUST_repo/main.nf \ 
             -resume \
             -profile conda \ 
             -c 16S.config \
             --reads 'sample.fastq'
```

See usage and output sections in the documentation (/docs) for all of the available options when running the pipeline.

## Computing requirements note

Clustering step uses up to 32-36GB RAM when working with a real dataset analysis and default parameters (umap_set_size = 100000). Setting umap_set_size to 50000, will diminish memory consumption to 10-13GB RAM (try using low_resource profile, eg. -profile conda,low_resource). When running the pipeline, kmer_freqs or mostly read_clustering processes could be terminated with status 137 when not enough RAM.

Nextflow automatically uses all available resources in your machine. More cpu threads enable the pipeline to compute and classify the different clusters at the same time and hence reduces the overall execution time.

Using the -with-trace option, it is possible to get an execution trace file which includes computing times and memory consumption metrics for all pipeline processes.

*The execution of the test profile (minimum testing dataset and default parameters) can be done with a regular 4 cores and 16GB RAM machine.

## Troubleshooting

- In some machines, the read_clustering process exits with error status(_RuntimeError: cannot cache function '...'_). We have seen that this condition can be avoided running the pipeline with sudo privileges (even if Docker was previously available without sudo permissions). 

## Credits

Rodríguez-Pérez H, Ciuffreda L, Flores C (2020). NanoCLUST: a species-level analysis of 16S rRNA nanopore sequencing data. Bioinformatics (2021) https://doi.org/10.1093/bioinformatics/btaa900

This work was supported by Instituto de Salud Carlos III [PI14/00844, PI17/00610, and FI18/00230] and co-financed by the European Regional Development Funds, “A way of making Europe” from the European Union; Ministerio de Ciencia e Innovación [RTC-2017-6471-1, AEI/FEDER, UE]; Cabildo Insular de Tenerife [CGIEU0000219140]; Fundación Canaria Instituto de Investigación Sanitaria de Canarias [PIFUN48/18]; and by the agreement with Instituto Tecnológico y de Energías Renovables (ITER) to strengthen scientific and technological education, training, research, development and innovation in Genomics, Personalized Medicine and Biotechnology [OA17/008]. 

## Contributions and Support

If you would like to contribute to this pipeline, please see the contributing guidelines
