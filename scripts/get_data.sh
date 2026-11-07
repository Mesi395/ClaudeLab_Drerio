#!/usr/bin/env bash
# Downloads the 22 public 10x libraries behind the problem dataset (Hu et al. 2022)
# into data/raw/. Source: GEO GSE159032 and GSE158919.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
raw="$root/data/raw"
mkdir -p "$raw"

geo="https://ftp.ncbi.nlm.nih.gov/geo/samples"

# GSM accession, GEO file label, library name used on the problem slide
libraries=(
  "GSM4817928 H5 H5"
  "GSM4817929 H6 H6"
  "GSM4817930 H7 H7"
  "GSM4817931 H8a H8a"
  "GSM4817932 H8v H8v"
  "GSM4815756 Hr11 Hr11"
  "GSM4817944 Hr12 Hr12"
  "GSM4817951 Hr24 Hr24"
  "GSM4817954 Hr27 Hr27"
  "GSM4817933 Hr1 Hr1"
  "GSM4817945 Hr13 Hr13"
  "GSM4817934 Hr2 Hr2a"
  "GSM4817935 Hr2_morecells Hr2b"
  "GSM4817938 Hr6a Hr6a"
  "GSM4817939 Hr6v Hr6v"
  "GSM4817940 Hr7a Hr7a"
  "GSM4817941 Hr7v Hr7v"
  "GSM4815755 Hr8 Hr8"
  "GSM4817942 Hr9 Hr9"
  "GSM4817947 Hr19 Hr19"
  "GSM4817936 Hr3 Hr3"
  "GSM4817937 Hr4 Hr4"
)

for entry in "${libraries[@]}"; do
  read -r gsm label lib <<<"$entry"
  out="$raw/${lib}_filtered_matrix.h5"
  [[ -s "$out" ]] && continue
  url="$geo/${gsm:0:7}nnn/$gsm/suppl/${gsm}_${label}_filtered_matrix.h5"
  echo "downloading $lib"
  curl -fsSL --retry 3 -o "$out" "$url"
done

echo "done: $(ls "$raw" | wc -l | tr -d ' ') libraries in data/raw"
