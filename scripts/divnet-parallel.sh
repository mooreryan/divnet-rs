#!/bin/bash

if [[ $# -ne 6 ]]; then
  echo >&2 "usage: divnet-parallel.sh <dn-fit-model> <dn-bootstrap> <path/to/outbase> <threads> <reps> <config.toml>"
  exit 1
fi

dn_fit_model="$1"
dn_bootstrap="$2"
base="$3"
threads="$4"
reps="$5"
config="$6"

dn_model_dat="$base".dat
dn_model_csv="$base".model.csv
final_outfile="$base".csv

DN_MODEL_DAT="$dn_model_dat" DN_MODEL_CSV="$dn_model_csv" "$dn_fit_model" "$config"

if [[ $? -ne 0 ]]; then
  echo >&2 "dn-fit-model failed"
  exit 1
fi

parallel --ungroup -j "$threads" \
  DN_MODEL_DAT="$dn_model_dat" DN_REPLICATE={} DN_SEED={} DN_BOOT_CSV="${base}".boot_{}.csv \
  "$dn_bootstrap" "$config" \
  ::: $(seq -s' ' 0 $(("$reps" - 1)))

if [[ $? -ne 0 ]]; then
  echo >&2 "dn-bootstrap failed"
  exit 1
fi

# Concatenate the outfiles
cat "$dn_model_csv" "${base}".boot_*.csv >"$final_outfile" &&
  rm "$dn_model_dat" "$dn_model_csv" "${base}".boot_*.csv
