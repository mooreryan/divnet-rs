library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Check arguments.
if (length(args) == 0) {
  stop("pass the input file.n", call. = FALSE)
}

csv <- args[[1]]
csv_out <- args[[2]]

csv <- read_csv(csv, comment = "#")

x <- csv[, 1:2]
y <- apply(csv[, 3:ncol(csv)], 2, signif)
csv <- cbind(x, y)

write_csv(csv, csv_out)
