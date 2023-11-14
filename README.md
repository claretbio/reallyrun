# ClaretBio's REALLY Library Processing Software

This software is for the basic informatic processing of sequencing data generated using ClaretBio's [REALLY](https://www.claretbio.com/really-technology) library prep kit with or without using unique molecular identifiers (UMIs).

# Installation

This software can be installed as a python package using the command `pip install reallyrun`

# Usage

The basic analysis can be run with `really runsamples` when running on standard libraries or `really runsamples --umi` when running on libraries with UMIs. The software takes in raw fastqs and trims adapters, aligns to a user-specified reference transcriptome, and marks duplicates. For UMI aware demltiplexing of REALLY libraries please use our SRSLYumi python package (more info at https://github.com/claretbio/SRSLYumi)

In order to run, this software requires an installation of conda. For speed, we recommend [mamba](https://mamba.readthedocs.io/en/latest/installation.html) which is best installed from [mambaforge](https://github.com/conda-forge/miniforge#mambaforge). If you prefer to use standard conda, installation instructions can be found [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links). 

Required Arguments

    --starindex : a path to the STAR index if one exists

    OR

    --reference: a path to the reference genome you wish to have converted to a STAR transcriptome for alignment

    --gtf: a path to the GTF of genes that corresponds to the reference you are aligning to
    
    --refflat: a path to the refFlat of genes that correspond to the reference you are aligning to

    --ribosomal: a path to the rRNA interval list that corresponds to the reference you are aligning to

    --libraries or --libfile: the library IDs you would like analzed in comma separated format or the path to a file with one ID per line, repsectively

Optional Arguments

    --fastqdir : a path to the directory containing the raw fastqs you wish to process (if not specified, defaults to current working directory)
    
    --resultsdir : a path to the directory you would like the output to be in (if not specified, defaults to current working directory)

    --indexdir: a path to the directory where you would like the STAR index to be created (if not specified, defalts to current working directory). Use one of --starindex or --indexdir, not both

Helpful information about arguments

The library IDs provided should match the beginning of the fastq files. For example, the library ID for the fastq files named `lib1_R1.fastq.gz` and `lib1_R2.fastq.gz` would be `lib1`. This can be provided directly on the command line with a comma separated list: `--libraries lib1,lib2` or as a file that lists one library ID per line: `--libfile libfile.txt`.

Details on GTF format can be found [here](https://useast.ensembl.org/info/website/upload/gff.html)

Picard's [BedToIntervalList tool](https://gatk.broadinstitute.org/hc/en-us/articles/360036716091-BedToIntervalList-Picard-) can be used to generate the ribosomal interval list. An example interval list can be found [here](https://gist.github.com/slowkow/b11c28796508f03cdf4b).

 For an example refFlat file for GRCh38, see `refFlat.txt.gz` at https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/.


Example Command
```
really runsamples --fastqdir /home/user/fastqfiles \
--resultsdir /home/user/amazing-results \
--reference /home/user/data/hg38.fa \
--gtf /home/user/data/hg38.gtf \
--refflat /home/user/data/hg38_refflat.txt \
--ribosomal /home/user/data/hg38_rrna.interval_list \
--libraries lib1,lib2,lib3
```

For reproducibility's sake and to ensure appropriate versions we use snakemake wrappers for many of the tools in this pipeline, which are often slow to create the first time they are used. As a result, your first time running the software may take a long time - don't worry, this is totally normal!
