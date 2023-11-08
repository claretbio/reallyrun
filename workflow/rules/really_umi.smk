rule fgbio_consensus:
    version:
        1
    message:
        "Consensus reads from UMIs with fgbio on {input}"
    output:
        "{base}/{sample}/{sample}.fgbio.cons.bam"
    input:
        "{base}/{sample}/{sample}.fgbio.bam"
    params:
        extra="-M 3 --tag=RX"
    log:
        "{base}/{sample}/{sample}-fgbio_consensus.log"
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
    wrapper:
        "v2.6.0/bio/fgbio/callmolecularconsensusreads"


rule fgbio_group:
    version:
        1
    output:
        bam="{base}/{sample}/{sample}.fgbio.bam",
        hist="{base}/{sample}/{sample}.fgbio.histo.tsv",
    input:
        "{base}/{sample}/{sample}.ppmate.bam"
    params:
        extra="--strategy=adjacency --min-map-q=0"
    log:
        "{base}/{sample}/{sample}-fgbio_group.log"
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
    wrapper:
        "v2.6.0/bio/fgbio/groupreadsbyumi"


rule filter_bam:
    version:
        1
    message:
        "Retaining only reads mapped in proper pair for {input}"
    output:
        bam="{base}/{sample}/{sample}.ppmate.bam",
    input:
        "{base}/{sample}/{sample}.rxmate.bam",
    params:
        extra="-f2",
    log:
        "{base}/{sample}/{sample}-samtools_view.log",
    threads: 4
    wrapper:
        "v2.6.0/bio/samtools/view"



rule fgbio_setmate:
    version:
        1
    message:
        "Setting mate tag with fgbio for {input}"
    output:
        "{base}/{sample}/{sample}.rxmate.bam", # was temp
    input:
        "{base}/{sample}/{sample}.rxsort.bam",
    params:
        extra=" -x",
    log:
        "{base}/{sample}/{sample}-fgbio_setmate.log",
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
    wrapper:
        "v2.6.0/bio/fgbio/setmateinformation"


rule query_name_sort:
    version:
        1
    message:
        "Sort {input} by query name",
    output:
        "{base}/{sample}/{sample}.rxsort.bam", # was temp
    input:
        "{base}/{sample}/{sample}.rx.dupbam",
    params:
        sort_order="queryname",
        extra="--VALIDATION_STRINGENCY LENIENT",
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
    wrapper:
        "v2.6.0/bio/picard/sortsam"


rule umi_tools_group:
    version:
        1
    message:
        "Correcting and grouping UMIs with UMI-Tools on {input}"
    output:
        "{base}/{sample}/{sample}.umi.dupbam",
        "{base}/{sample}/{sample}.umitools.tsv",
    input:
        "{base}/{sample}/{sample}.rx.dupbam",
        "{base}/{sample}/{sample}.rx.dupbam.bai",
    log:
        "{base}/{sample}/{sample}-umitools.log",
        "{base}/{sample}/{sample}-umitools.err"
    params:
        min_map_q=0,
        random_seed=69,
    shell:
        """
        umi_tools group --output-bam --stdin={input} --stdout={output[0]} \
        --group-out={output[1]} --log={log[0]} --error={log[1]} \
        --extract-umi-method=tag --umi-tag=RX --method=directional \
        --paired --unmapped-reads=use --mapping-quality={params.min_map_q} \
        --random-seed={params.random_seed}
        """

rule umi_srslyumi_bamtag:
    version:
        1
    message:
        "Moving UMI from fragment name to RX tag on {input}"
    output:
        "{base}/{sample}/{sample}.rx.dupbam", # was temp
    input:
        "{base}/{sample}/{sample}.dupbam",
        "{base}/{sample}/{sample}.dupbam.bai",
    log:
        "{base}/{sample}/{sample}-rxdupbam.log",
    shell:
        """
        srslyumi-bamtag --binary -o {output} --take-fragment 0 {input[0]} 2> {log}
        """
