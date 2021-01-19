LEE_PHYLUM_DIR="./test_files/lee_phylum"

# Try to delete output files.
rm "${LEE_PHYLUM_DIR}/lee_phylum_divnet_output.csv" \
  "${LEE_PHYLUM_DIR}/lee_phylum_divnet_plot.png"

"${DIVNET_RS}" "${LEE_PHYLUM_DIR}/config_lee_phylum.toml" && \
  "${RSCRIPT}" --vanilla "${LEE_PHYLUM_DIR}/scripts/make_plot.R" "${LEE_PHYLUM_DIR}"
