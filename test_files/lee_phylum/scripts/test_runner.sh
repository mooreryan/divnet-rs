LEE_PHYLUM_DIR="./test_files/lee_phylum"

"${DIVNET_RS}" "${LEE_PHYLUM_DIR}/config_lee_phylum.toml" && \
  "${RSCRIPT}" --vanilla "${LEE_PHYLUM_DIR}/scripts/make_plot.R" "${LEE_PHYLUM_DIR}"
