LEE_PHYLUM_DIR = "./test_files/lee_phylum"

.PHONY: test_lee_phylum

test_lee_phylum:
	${SH} $(LEE_PHYLUM_DIR)/scripts/test_runner.sh
