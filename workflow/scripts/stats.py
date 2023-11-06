#!/usr/bin/env python3

import argparse
import sys
import subprocess


ls = {
    "Number of input reads": "input_reads",
    "Average input read length": "avg_input_length",
    "Average mapped length": "avg_length_mapped",
    "Uniquely mapped reads number": "num_mapped_unique",
    "Uniquely mapped reads %": "perc_mapped_unique",
    "Number of splices: Total": "num_splice_aligned",
    "Number of reads mapped to multiple loci": "num_mapped_multi",
    "% of reads mapped to multiple loci": "perc_mapped_multi",
    "% of reads unmapped: too short": "perc_unmapped_too_short",
    "Number of chimeric reads": "num_chimeric",
}


def read_starlog(fn, ls):
    r = {}
    with open(fn) as f:
        for line in f:
            if "|" in line:
                sp = line.split("|")
                sp[0] = sp[0].strip(" ")
                sp[1] = sp[1].strip("\t,\n")
                if sp[0] in ls:
                    k = ls[sp[0]]
                    if "%" in sp[1]:
                        r[k] = sp[1].strip("%")
                    else:
                        r[k] = sp[1]

    return r


def read_flagstats(fn):
    r = {}
    with open(fn) as f:
        for line in f:
            count, _, _, *description = line.rstrip().split(" ")
            key = " ".join(description)
            if key.startswith("properly paired"):
                key = "properly paired"
            elif key.startswith("mapped"):
                key = "mapped"
            elif key.startswith("singletons"):
                key = "singletons"
            r[key] = int(count)
    r["total"] = (
        r["in total (QC-passed reads + QC-failed reads)"]
        - r["secondary"]
        - r["supplementary"]
    )
    r["fraction_pairs_mapped"] = (
        r["with itself and mate mapped"] + 2 * r["singletons"]
    ) / float(r["total"])
    r["fraction_chimeras"] = r["with mate mapped to a different chr"] / float(
        r["total"]
    )
    return r


def read_cutadapt(fn):
    r = {}
    with open(fn) as f:
        for line in f:
            if line.startswith("Total read pairs"):
                key, _, val = line.rstrip().partition(":")
                r[key] = val.strip()
            if line.startswith("Pairs"):
                key, _, val = line.rstrip().partition(":")
                r[key] = val.strip()
    return r


def read_picard_metrics_file(fn):
    conv = {"LIBRARY": str, "PERCENT_DUPLICATION": float}
    with open(fn) as f:
        for line in f:
            if line.startswith("## METRICS CLASS"):
                break
        header = f.readline().rstrip().split("\t")
        values = f.readline().rstrip().split("\t")
        return dict([(h, conv.get(h, float)(v)) for h, v in zip(header, values)])


def main(
    libname,
    cutadaptfn,
    starfn,
    picard_RNAmetrics,
    picard_dupMetrics,
    out,
):
    cutadapt = read_cutadapt(cutadaptfn)
    star = read_starlog(starfn, ls)
    picardRNA = read_picard_metrics_file(picard_RNAmetrics)
    picardDup = read_picard_metrics_file(picard_dupMetrics)
    o = [
        ("lib", libname),
        ("pairs", cutadapt["Total read pairs processed"]),
        ("pairs_too_short", cutadapt["Pairs that were too short"]),
        ("pairs_kept", cutadapt["Pairs written (passing filters)"]),
        ("avg_input_length", star["avg_input_length"]),
        ("avg_length_mapped", star["avg_length_mapped"]),
        ("mapped_unique", star["num_mapped_unique"]),
        ("percent_mapped_unique", float(star["perc_mapped_unique"])),
        ("num_splice_aligned", star["num_splice_aligned"]),
        ("num_multimappers", star["num_mapped_multi"]),
        ("percent_multimappers", float(star["perc_mapped_multi"])),
        ("percent_unmapped_tooshort", float(star["perc_unmapped_too_short"])),
        ("num_chimeric", star["num_chimeric"]),
        (
            "picard_percent_dup",
            "{:.2f}".format((picardDup["PERCENT_DUPLICATION"]) * 100),
        ),
        ("picard_est_lib_size", (picardDup["ESTIMATED_LIBRARY_SIZE"])),
        ("ribosomal_bases", picardRNA["RIBOSOMAL_BASES"]),
        (
            "percent_ribosomal_bases",
            "{:.2f}".format(picardRNA["PCT_RIBOSOMAL_BASES"] * 100),
        ),
        ("correct_strand", picardRNA["CORRECT_STRAND_READS"]),
        (
            "percent_correct_strand",
            "{:.2f}".format(picardRNA["PCT_CORRECT_STRAND_READS"] * 100),
        ),
        ("incorrect_strand", picardRNA["INCORRECT_STRAND_READS"]),
        ("percent_mrna_bases", "{:.2f}".format(picardRNA["PCT_MRNA_BASES"] * 100)),
        ("percent_coding_bases", "{:.2f}".format(picardRNA["PCT_CODING_BASES"] * 100)),
        ("percent_utr_bases", "{:.2f}".format(picardRNA["PCT_UTR_BASES"] * 100)),
        (
            "percent_intronic_bases",
            "{:.2f}".format(picardRNA["PCT_INTRONIC_BASES"] * 100),
        ),
        (
            "percent_intergenic_bases",
            "{:.2f}".format(picardRNA["PCT_INTERGENIC_BASES"] * 100),
        ),
        ("percent_usable_bases", "{:.2f}".format(picardRNA["PCT_USABLE_BASES"] * 100)),
    ]

    for key, val in o:
        out.write("{}\t{}\n".format(key, val))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "-o",
        default=sys.stdout,
        help="output file (stdout)",
        type=argparse.FileType("w"),
    )
    ap.add_argument("libname", help="Library name and seq date")
    ap.add_argument("cutadapt", help="cutadapt output")
    ap.add_argument("star", help="STAR final log output")
    ap.add_argument("picard_rnametrics", help="Picard tools RNA metrics")
    ap.add_argument("picard_dupmetrics", help="Picard tools duplicate metrics")

    a = ap.parse_args()
    main(
        a.libname,
        a.cutadapt,
        a.star,
        a.picard_rnametrics,
        a.picard_dupmetrics,
        a.o,
    )
