import glob

rule stats:
    version:
        1
    message:
        "Collating stats to {output}"
    output:
        "{resultsdir}/{{sample}}/{{sample}}.stats".format(**config),
    input:
        cutadapt="{resultsdir}/{{sample}}/{{sample}}.qc.txt".format(**config),
        starlog="{resultsdir}/{{sample}}/{{sample}}.Log.final.out".format(**config),
        flagstat="{resultsdir}/{{sample}}/{{sample}}.bam.flagstat".format(**config),
        picardRNA="{resultsdir}/{{sample}}/{{sample}}.RNAseq_metrics.txt".format(**config),
        picardDupe="{resultsdir}/{{sample}}/{{sample}}.markdupe_metrics.txt".format(**config),
    shell:
        """
        python3 workflow/scripts/stats.py {wildcards.sample} {input.cutadapt} {input.starlog} \
        {input.picardRNA} {input.picardDupe} -o {output}
        """


rule gcbias:
    version:
        1
    message:
        "Calculating and plotting GC bias for {input}"
    output:
        xls="{base}/{sample}/{sample}.GC.xls",
        r=temp("{base}/{sample}/{sample}.GC_plot.r"),
        pdf="{base}/{sample}/{sample}.GC_plot.pdf",
    input:
        "{base}/{sample}/{sample}.bam",
    params:
        outpre="{base}/{sample}/{sample}",
    conda:
        "../envs/r-script.yml"
    shell:
        """
        python3 workflow/scripts/read_GC.py -i {input} -o {params.outpre}
        """

rule flagstat:
    version:
       1
    message:
       "Running samtools flagstat on {input}"
    output:
        "{base}/{sample}/{sample}.bam.flagstat",
    input:
        "{base}/{sample}/{sample}.bam",
    log: 
        "{base}/{sample}/{sample}.flagstat.log"
    wrapper:
        "v2.6.0/bio/samtools/flagstat"


rule plot_mol_len:
    version:
        1
    output:
        "{base}/{sample}/{sample}.insert_length_dist.pdf"
    input:
        "{base}/{sample}/{sample}.mol_lens.txt"
    conda:
        "../envs/r-script.yml"
    shell:
        "Rscript workflow/scripts/plot_mol_len.R {input} $(basename {wildcards.base}) {output}" 


rule write_mol_len:
    version:
        1
    output:
        temp("{base}/{sample}/{sample}.mol_lens.txt")
    input:
        bam="{base}/{sample}/{sample}.mml.bam",
    shell:
        "python3 workflow/scripts/write_mol_len.py {input.bam} {output}"


rule filter_mol_len:
    version:
        1
    output:
        bam=temp("{base}/{sample}/{sample}.mml.bam"),
    input:
        "{base}/{sample}/{sample}.bam"
    params:
        extra="-f64 -F2308",
        samtools_opts="--write-index"
    wrapper:
        "v2.6.0/bio/samtools/view"


rule rna_seq_metrics:
    version:
        1
    message:
        "Running picard on {input} to get RNA seq metrics",
    output:
        "{base}/{sample}/{sample}.RNAseq_metrics.txt",
    input:
        bam="{base}/{sample}/{sample}.bam",
        refflat=config.get("refflat")
    params:
        strand="SECOND_READ_TRANSCRIPTION_STRAND",
        extra="--VERBOSITY ERROR --RIBOSOMAL_INTERVALS {}".format(config.get("ribosomal"))
    log:
        "{base}/{sample}/{sample}.RNAseq_metrics.log",
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
    wrapper:
        "v2.6.0/bio/picard/collectrnaseqmetrics"


rule bam_index:
    version:
        1
    message:
        "Indexing bam {input}"
    output:
        "{base}bam.bai",
    input:
        "{base}bam",
    log:
        "{base}.samtools_index.log",
    wrapper:
        "v2.6.0/bio/samtools/index"


def get_params(umi):
    if umi == 'true':
        return "--BARCODE_TAG BX --REMOVE_DUPLICATES true"
    else:
        return "--REMOVE_DUPLICATES true"


def get_mark_dup_input(umi):
    if umi == 'true':
        return "{base}/{sample}/{sample}.umi.dupbam"
    else:
        return "{base}/{sample}/{sample}.dupbam"


rule mark_dupe:
    version:
        "1"
    message:
        "Marking duplicates using Picard",
    output:
        bam="{base}/{sample}/{sample}.bam",
        metrics="{base}/{sample}/{sample}.markdupe_metrics.txt",
    input:
        bams=get_mark_dup_input(config['umi'])
    log:
        "{base}/{sample}/{sample}.markdupe.log",
    params:
        extra=get_params(config['umi'])
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
    wrapper:
        "v2.6.0/bio/picard/markduplicates"


def get_star_index():
    if config.get("starindex") is None:
        return "{}/starIndex".format(config.get("indexdir"))
    else:
        return config.get("starindex")


rule star_align:
    version:
        "1"    
    message:
        "Aligning trimmed fastqs with STAR",
    output:
        aln="{base}/{sample}/{sample}.dupbam",
        log_final="{base}/{sample}/{sample}.Log.final.out",
        reads_per_gene="{base}/{sample}/{sample}.ReadsPerGene.out.tab",
        sj="{base}/{sample}/{sample}.SJ.out.tab",
        unmapped=["{base}/{sample}/{sample}_R1.unmapped.fq.gz",
                    "{base}/{sample}/{sample}_R2.unmapped.fq.gz"],
    input:
        fq1="{base}/{sample}/{sample}_R1.fastq",
        fq2="{base}/{sample}/{sample}_R2.fastq",
        idx=get_star_index(),
    log:
        "{base}/{sample}/{sample}.star.log",
    params:
        extra="--sjdbGTFfile {gtf} --outSAMtype BAM SortedByCoordinate  --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 --outFilterMismatchNmax 2 --bamRemoveDuplicatesType UniqueIdenticalNotMulti --quantMode GeneCounts".format(**config)
    threads: 8
    wrapper:
        "v2.6.0/bio/star/align"


class TooManyMatchingFastqError(Exception):
    def __init__(self, read, base, rawdir, files):
        self.read = read
        self.base = base
        self.rawdir = rawdir
        self.files = files


class FastqMissingError(Exception):
    def __init__(self, read, base, rawdir):
        self.read = read
        self.base = base
        self.rawdir = rawdir


def find_a_fastq(indir, sample, read):
    pattern = "{}/{}*{}*.fastq.gz".format(indir, sample, read)
    reads = glob.glob(pattern)

    if len(reads) == 0:
        raise FastqMissingError(indir, sample, read)
    elif len(reads) > 1:
        raise TooManyMatchingFastqError(indir, sample, read, reads)

    return reads[0]


def find_fastqs(sample):
    r1 = find_a_fastq(config["indir"], sample, "R1")
    r2 = find_a_fastq(config["indir"], sample, "R2")
    return [r1, r2]


rule cutadapt:
    version:
        1    
    message:
        "Trimming adapters with cutadapt"
    output:
        fastq1="{resultsdir}/{{sample}}/{{sample}}_R1.fastq".format(**config),
        fastq2="{resultsdir}/{{sample}}/{{sample}}_R2.fastq".format(**config),
        qc="{resultsdir}/{{sample}}/{{sample}}.qc.txt".format(**config)
    input:
        lambda wildcards: find_fastqs(wildcards.sample)
    params:
        adapters="-a AGATCGGAAGAGCACACGTCTGAA -g AGATCGGAAGAGCGTCGTGTAGGG -A AGATCGGAAGAGCACACGTCTGAA -G AGATCGGAAGAGCGTCGTGTAGGG",
        extra="--minimum-length 30"
    log:
       "{resultsdir}/{{sample}}/{{sample}}-cutadapt.log".format(**config)
    wrapper:
        "v2.6.0/bio/cutadapt/pe"


rule star_index:
    version:
        1
    input:
        fasta="{}".format(config.get("reference")),
    output:
        directory("{}/starIndex".format(config.get("indexdir"))),
    message:
        "Creating STAR index"
    threads: 8
    params:
        extra="--sjdbGTFfile {}".format(config.get("gtf")), ## add gtf file
    log:
        "{}/starIndex/star_index.log".format(config.get("indexdir")),
    wrapper:
        "v2.6.0/bio/star/index"
