library(tidyverse)
library(DivNet)
library(tictoc)

data(Lee)

set.seed(1537290)

# You may need to change the dir.
#
# It is a relative path from where you're running the script!
dir <- "./docs/lee_data"
plot_out <- file.path(dir, "lee_phylum_divnet_plot.png")
count_table_out <- file.path(dir, "lee_phylum_count_table.csv")
sample_data_out <- file.path(dir, "lee_phylum_sample_data.csv")

# Collapse taxa by phyla.
lee <- tax_glom(Lee, taxrank = "Phylum")

# Run DivNet.
tic("divnet lee phylum")
dn <- divnet(lee,
             X = "char",
             ncores = 1,
             network = "diagonal",
             base = 1,
             B = 5,
             tuning = "fast")
toc()

# Plot results.
png(plot_out,
    600,
    600,
    pointsize = 24)
plot(dn) + ggtitle("DivNet")
dev.off()

# Make sure taxa are rows.
stopifnot(taxa_are_rows(Lee))

# Write the count table.
lee %>%
  otu_table %>%
  as.data.frame %>%
  rownames_to_column("taxa") %>%
  write.table(count_table_out,
              quote = FALSE,
              sep = ",",
              row.names = FALSE)

# Write the sample data.
#
# There is a bit of hacking you have to do....
d <- lee %>% sample_data

class(d) <- "data.frame"

model.matrix(~char, d)[, 2:5] %>%
  as.data.frame %>%
  rownames_to_column("sample") %>%
  write.table(sample_data_out,
              quote = FALSE,
              sep = ",",
              row.names = FALSE)
