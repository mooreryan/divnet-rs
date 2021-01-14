library(tidyverse)
library(DivNet)

data(Lee)

# You will need to change these!
count_table_fname <- "lee_count_table.csv"
sample_data_fname <- "lee_sample_data.csv"

if (taxa_are_rows(Lee)) {
  dat <- Lee %>%
    otu_table %>%
    as.data.frame %>%
    rownames_to_column("taxa")

  dat %>%
    write.table(count_table_fname,
                quote = FALSE,
                sep = ",",
                row.names = FALSE)
} else {
  # TODO just transpose it!
  stop("whoops, taxa aren't rows")
}


d <- Lee %>% sample_data

class(d) <- "data.frame"

model.matrix(~char, d)[, 2:5] %>%
  as.data.frame %>%
  rownames_to_column("sample") %>%
  write.table(sample_data_fname,
              quote = FALSE,
              sep = ",",
              row.names = FALSE)
