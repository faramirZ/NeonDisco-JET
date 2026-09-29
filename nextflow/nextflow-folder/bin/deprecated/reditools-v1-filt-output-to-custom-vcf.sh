#!/usr/bin/env bash
# reditools-v1-filt-output-to-custom-vcf.sh
# Convert REDItools2 output TSV to VCF format for use with SnpEff
#
# Usage: ./reditools-v1-filt-output-to-custom-vcf.sh input.tsv output.vcf

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input.tsv output.vcf" >&2
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

if [[ ! -f "$INPUT" ]]; then
    echo "Error: input file '$INPUT' not found" >&2
    exit 1
fi

awk -F'\t' 'BEGIN{
    print "##fileformat=VCFv4.1"
    print "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand\">"
    print "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"RNA coverage\">"
    print "##INFO=<ID=MEANQ,Number=1,Type=Float,Description=\"Mean quality\">"
    print "##INFO=<ID=BASECOUNT,Number=1,Type=String,Description=\"Base counts [A,C,G,T]\">"
    print "##INFO=<ID=FREQ,Number=1,Type=Float,Description=\"Editing frequency\">"
    print "##INFO=<ID=GCOV,Number=1,Type=Integer,Description=\"Genomic coverage\">"
    print "##INFO=<ID=GMEANQ,Number=1,Type=Float,Description=\"Genomic mean quality\">"
    print "##INFO=<ID=GBASECOUNT,Number=1,Type=String,Description=\"Genomic base counts [A,C,G,T]\">"
    print "##INFO=<ID=GALLSUBS,Number=1,Type=String,Description=\"Genomic substitutions\">"
    print "##INFO=<ID=GFREQ,Number=1,Type=Float,Description=\"Genomic frequency\">"
    print "##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Repeat type\">"
    print "##INFO=<ID=REPNAME,Number=1,Type=String,Description=\"Repeat name\">"
    print "##INFO=<ID=SNP,Number=1,Type=String,Description=\"SNP flag\">"
    print "##INFO=<ID=DBSNP,Number=1,Type=String,Description=\"dbSNP ID\">"
    print "##INFO=<ID=REDIPORTAL,Number=1,Type=String,Description=\"REDIPortal known editing site\">"
    print "##INFO=<ID=EDITSTATUS,Number=1,Type=String,Description=\"Editing status\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

    # Complement lookup for minus strand correction
    comp["A"] = "T"; comp["T"] = "A"
    comp["C"] = "G"; comp["G"] = "C"
}
NR==1{next}
{
    ref = $3
    alt = substr($8, 2, 1)
    if(alt == "." || alt == "") next

    # For minus strand sites (strand=0 in REDItools v1), flip REF and ALT to genomic forward strand
    if ($4 == "0") {
        ref = comp[ref]
        alt = comp[alt]
    }

    strand   = ($4  == "-" ? "." : $4)
    reptype  = ($15 == "-" ? "." : $15)
    repname  = ($16 == "-" ? "." : $16)
    snpflag  = ($17 == "-" ? "." : $17)
    dbsnp    = ($18 == "-" ? "." : $18)
    redi     = ($19 == "-" ? "." : $19)
    editstat = ($20 == "-" ? "." : $20)
    gallsubs = ($13 == "-" ? "." : $13)

    info = sprintf("STRAND=%s;COV=%s;MEANQ=%s;BASECOUNT=%s;FREQ=%s;GCOV=%s;GMEANQ=%s;GBASECOUNT=%s;GALLSUBS=%s;GFREQ=%s;REPTYPE=%s;REPNAME=%s;SNP=%s;DBSNP=%s;REDIPORTAL=%s;EDITSTATUS=%s",
        strand, $5, $6, $7, $9, $10, $11, $12, gallsubs, $14, reptype, repname, snpflag, dbsnp, redi, editstat)

    print $1"\t"$2"\t.\t"ref"\t"alt"\t.\tPASS\t"info
}' OFS='\t' "$INPUT" > "$OUTPUT"

echo "Done: $OUTPUT"
