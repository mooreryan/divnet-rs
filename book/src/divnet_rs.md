# divnet-rs

[DivNet](https://github.com/adw96/DivNet) is an R package for estimating diversity when taxa in the community co-occur via ecological networks.  It leverages information from multiple samples and data about the samples (covariates) in order to give more accurate estimates of variance in the measured diversity.  Possibly it's coolest feature is that, unlike most existing methods of measuring diversity, it uses models from [compositional data analysis](https://en.wikipedia.org/wiki/Compositional_data) that take into account ecological networks in which taxa positively and negatively co-occur.

While the DivNet R package makes it simple to apply the algorithms from [Willis & Martin 2020](https://github.com/adw96/DivNet) to your data, it does have trouble with large datasets.  This is where `divnet-rs` comes it.  It is a [Rust](TODO) implementation of the DivNet algorithm.  It allows you to successfully run much larger datasets even on your laptop!

In this book, you will find documentation and information for installing and using `divnet-rs`.  Additionally, I will go over a couple examples and a tutorial to help you get started!

For more background information and the theory behind DivNet, please check out the [original DivNet repository](https://github.com/adw96/DivNet) and the [DivNet manuscript](https://doi.org/10.1093/biostatistics/kxaa015).

You can find the `divnet-rs` source code on [GitHub](https://github.com/mooreryan/divnet-rs).
