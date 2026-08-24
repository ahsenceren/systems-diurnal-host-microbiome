#!/usr/bin/env bash
# 01_fetch_metagenome.sh
set -euo pipefail

PROJECT="PRJEB89098"
OUTDIR="data/raw/metagenome"
mkdir -p "$OUTDIR"

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,fastq_ftp,library_strategy,library_source,sample_title" \
    > "${OUTDIR}/all_runs_manifest.tsv"

echo "Toplam run sayısı:"
wc -l < "${OUTDIR}/all_runs_manifest.tsv"

echo ""
echo "library_strategy dağılımı (metaT vs metaG ayrımını görmek için):"
cut -f3 "${OUTDIR}/all_runs_manifest.tsv" | sort | uniq -c
