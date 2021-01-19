library(tictoc)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

# Check arguments.
if (length(args) == 0) {
  stop("pass the input dir in.n", call. = FALSE)
}

dir <- args[[1]]
plot_out <- file.path(dir, "lee_phylum_divnet_plot.png")
divnet_rs_in <- file.path(dir, "lee_phylum_divnet_output.csv")

# Sorry y'all, currently, you still have to set this manually!
nreplicates <- 5

# Read data file.
divnet_rs <- read.table(divnet_rs_in, sep = ",", header = TRUE)

# Replicate 0 is actually the estimates for the real data.
rep0 <- divnet_rs[divnet_rs$replicate == 0, -1]
rownames(rep0) <- rep0$sample
rep0$sample <- NULL

# Get the Shannon index for the actual data.
rep0_shannon <- apply(rep0, 1, DivNet::shannon_true)

# Now calculate the shannon index for the replicates.
reps_shannon <- sapply(1:nreplicates, function (i) {
  d <- divnet_rs[divnet_rs$replicate == i, -1]
  rownames(d) <- d$sample
  d$sample <- NULL

  apply(d, 1, DivNet::shannon_true)
})

# What we want is the variance in the diversity estimates for the replicates.
reps_shannon_error <- t(apply(reps_shannon, 1, function (x) {
  c(var(x), sd(x))
}))
colnames(reps_shannon_error) <- c("variance", "sd")


# Make the plot!  To be fair, this could use some cleaning up!
png(plot_out,
    600,
    600,
    pointsize = 24)
tibble(names = names(rep0_shannon),
       shannon = rep0_shannon) %>%
  left_join(reps_shannon_error %>%
              as.data.frame %>%
              rownames_to_column("names")) %>%
  mutate(ciupper = shannon + 2 * sd,
         cilower = shannon - 2 * sd) %>%
  ggplot(aes(x = names)) +
  geom_segment(aes(xend = names, y = cilower, yend = ciupper)) +
  geom_point(aes(y = shannon)) +
  xlab("samples") +
  ylab("shannon estimate") +
  ggtitle("divnet-rs") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.off()
