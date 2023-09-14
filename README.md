# wf-emu

This repository contains a [nextflow](https://www.nextflow.io/) workflow for creating a microbial community profile from 16s rRNA targeted amplicon sequences. 


## Introduction

`wf-emu` utilize [Emu](https://gitlab.com/treangenlab/emu) for creation of microbial community profile. Emu uses [minimap2](https://github.com/lh3/minimap2) to align sample reads to a database, then employs an expectation-maximization algorithm to reweight read classifications based on relative abundance estimates within the sample. Reads are expected to be demultiplexed (of a single sample) and barcodes are expected to be removed prior to running this workflow. Additional quality control is up to the user's discretion.


## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run wf-emu --help
```

to see the options for the workflow.

The main options for the workflow are:
* `fastq`: A fastq file or directory containing fastq input files or directories of input files.
* `analyse_unclassified`: Set to true to generate additional output file of unclassified sequences.
* `sample_sheet`: A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. 
* `sample`: A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files.
* `out_dir`: Specify directory for all workflow results.
* `min_abundance`: Minimum species relative abundance in results. Species with estimated abundance below this value will be marked as not present.
* `K`: Minibatch size for minimap2 mapping.
* `threads`: Number of CPU threads during minimap2 alignment step.
* `disable_ping`: enable to prevent sending a workflow ping.


**Workflow outputs**

The primary outputs of the workflow include:

* rel_abun directory: .tsv file with relative abundance for each sample.
* read_assignment directory: .tsv file for each sample that contains classification distributions for each read within the sample.
* unclassified directory: fasta file with unclassified sequences for each sample.
* emu-combined-x.tsv: single file to combine community profiles for each sample included in the workflow.
* wf-emu-report.html: page with all read quality information and community profiles (with sunburst plot visual) for each sample in the workflow.

**Databases**

`wf-emu` is set up with one default database. This is a summation of [NCBI 16S RefSeq](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/) from September 2020 and [rrnDB v5.6](https://rrndb.umms.med.umich.edu/). Taxonomy is from NCBI on the same download date. The resulting database contains 49,301 sequences from 17,555 unique bacterial and archaeal species.



## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
