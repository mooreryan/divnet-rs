LEE_PHYLUM_DIR="./test_files/lee_phylum"
PREVIOUS_OUTPUT="${LEE_PHYLUM_DIR}/expected_output"

# Try to delete output files.
rm "${LEE_PHYLUM_DIR}/lee_phylum_divnet_output.csv" \
  "${LEE_PHYLUM_DIR}/lee_phylum_divnet_plot.png"

"${DIVNET_RS}" "${LEE_PHYLUM_DIR}/config_lee_phylum.toml" && \
  "${RSCRIPT}" --vanilla "${LEE_PHYLUM_DIR}/scripts/make_plot.R" "${LEE_PHYLUM_DIR}" && \
  "${RSCRIPT}" --vanilla "${LEE_PHYLUM_DIR}/scripts/round.R" \
   "${PREVIOUS_OUTPUT}/lee_phylum_divnet_output.csv" "${PREVIOUS_OUTPUT}/lee_phylum_divnet_output.csv.ROUND" && \
 "${RSCRIPT}" --vanilla "${LEE_PHYLUM_DIR}/scripts/round.R" \
   "${LEE_PHYLUM_DIR}/lee_phylum_divnet_output.csv" "${LEE_PHYLUM_DIR}/lee_phylum_divnet_output.csv.ROUND" && \
  diff "${PREVIOUS_OUTPUT}/lee_phylum_divnet_output.csv.ROUND" "${LEE_PHYLUM_DIR}/lee_phylum_divnet_output.csv.ROUND" && \
  diff "${PREVIOUS_OUTPUT}/lee_phylum_divnet_plot.png" "${LEE_PHYLUM_DIR}/lee_phylum_divnet_plot.png"
