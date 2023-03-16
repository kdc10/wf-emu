# Emu workflow for EPI2ME

## Introduction

This repository contains a [nextflow](https://www.nextflow.io/) derived from 
[template](https://github.com/epi2me-labs/wf-template) for [Emu](https://gitlab.com/treangenlab/emu).
Emu is a community profile estimator for 16S rRNA amplicon sequences. 
The method is optimized for error-prone full-length reads, but can also be utilized for short-read data.

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

To test installation, users can run:

```
nextflow run wf-emu --fastq test_data/Zymo_ONT.fastq
```

**Workflow outputs**

The primary outputs of the workflow include:

* a simple text file providing a summary of sequencing reads,
* an HTML report document detailing the primary findings of the workflow.
* Emu-derived community profiles 

## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
