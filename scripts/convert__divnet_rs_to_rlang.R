library(DivNet)
library(tictoc)
library(tidyverse)

# Of course, you will probably want to handle the sample data as well!

# You will need to change these!
dn_rs_infile <- "./test_files/lee_phylum/lee_phylum_divnet_output.csv"
nreplicates <- 5 # This should match your divnet-rs config!
png_output <- "lee_divnet_shannon.png"

dn_rs <- read.table(dn_rs_infile, sep = ',', header = T)

dn_rs_rep0 <- dn_rs[dn_rs$replicate == 0, -1]
rownames(dn_rs_rep0) <- dn_rs_rep0$sample
dn_rs_rep0$sample <- NULL

# Get the shannon variance.....each row is a sample, each col is the taxa.
shan_replicates <- sapply(1:nreplicates, function(i) {
  d <- dn_rs[dn_rs$replicate == i, -1]
  rownames(d) <- d$sample
  d$sample <- NULL

  shan <- apply(d, 1, DivNet:::shannon_true)
})

shan_error <- t(apply(shan_replicates, 1, function(x) {
  c(var(x), sd(x))
}))
colnames(shan_error) <- c("variance", "sd")

shan_error <- shan_error %>%
  as_tibble %>%
  mutate(names = rownames(shan_replicates))

shannon <- apply(dn_rs_rep0, 1, DivNet:::shannon_true)

plotting <- tibble(x = 1:length(shannon), shannon = shannon, names = names(shannon))

size <- 600
pointsize <- 24

png(png_output, size, size, pointsize = pointsize)
plotting %>%
  left_join(shan_error) %>%
  mutate(ciupper = shannon + 2 * sd,
         cilower = shannon - 2 * sd) %>%
  ggplot(aes(x = x)) +
  scale_x_continuous(breaks = plotting$x, labels = plotting$names) +
  geom_segment(aes(x = x, xend = x, y = cilower, yend = ciupper)) +
  geom_point(aes(y = shannon)) +
  xlab("samples") +
  ylab("shannon estimate") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.off()
