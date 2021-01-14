# divnet-rs

It's [DivNet](https://github.com/adw96/DivNet) with fewer features, but also *way* faster!!

Unfortunately, it is also more complicated to use at the moment, so if your data sets are smaller, I recommend using the [original implementation](https://github.com/adw96/DivNet).

Please see the [GitHub repository for the reference implementation](https://github.com/adw96/DivNet).  It has lots of good info about the theory behind DivNet!

## Installation

### Precompiled binaries

TODO

### From source

`divnet-rs` is written in [Rust](https://www.rust-lang.org/).  If you do not have Rust, you must [install it first](https://www.rust-lang.org/tools/install).

There are two options.  You can either follow the master branch (recommended only if you want to do development) or download a tagged release.

If you want to track the master branch, clone the git repository like so:

```
git clone https://github.com/mooreryan/divnet-rs.git
```

I recommend you to download a [release](TODO) however.

```
TODO
```

Once you have the source code, `cd` in to the root of the source directoy and run `cargo` to compile the program.  Here is an example:

```
cd divnet-rs
cargo build --release
```

Note that this can take a while, but you will *not* have to compile the program each time you use it!

You should now have a compiled binary in `./target/release/divnet-rs`.  You will probably want move or symlink the executable binary (`./target/release/divnet-rs`) somewhere on your path.  For example, if you want to create a symlink to `~/bin`, you could run something like this:

```
ln -s $(pwd)/target/release/divnet-rs $HOME/bin/divnet-rs
```

You can check that you have access to the `divnet-rs` binary with the `which` command:

```
which divnet-rs #=> /home/ryan/bin/divnet-rs
```

If you see a path, then the linking worked!  If not, then you can still just use the program from its current directory (`target/release/divnet-rs`) or try and [figure out how to get it on your path](https://gist.github.com/nex3/c395b2f8fd4b02068be37c961301caa7).

### Testing your installation

The `divnet-rs` repository includes some handy test files to check that everything is working correctly.

Try it out!

```
divnet-rs ./test_files/small/config_small.toml
```

If all goes well, you should see some logging output that ends with a line like this:

```
[2021-01-14T02:59:49Z INFO  divnet_rs] Done!  Output file is: "./test_files/small/small_divnet_output.csv"
```

## Usage

See the `./docs` (TODO) directory included in this repository for information on how to run `divnet-rs`, how to get your data from R (and possibly phyloseq) into a format that `divnet-rs` can handle, and how to get your data *back* into R so you can get on with your research!

## Citation

If you use DivNet, please support the authors and [cite their manuscript](https://doi.org/10.1093/biostatistics/kxaa015!

## Important differences from reference implementation

### Diagonal network only

The only network option is `diagonal`.  Here is an [quote from Amy Willis](https://github.com/adw96/DivNet/issues/32#issuecomment-521727997) that supports this decision:

> I would recommend network="diagonal" for a dataset of this size. This means you're allowing overdispersion (compared to a plugin aka multinomial model) but not a network structure. This isn't just about computational expense -- it's about the reliability of the network estimates. Essentially estimating network structure on 20k variables (taxa) with 50 samples with any kind of reliability is going to be very challenging, and I don't think that it's worth doing here. In our simulations we basically found that overdispersion contributes the bulk of the variance to diversity estimation (i.e. overdispersion is more important than network structure), so I don't think you are going to lose too much anyway.

The [whole issue](https://github.com/adw96/DivNet/issues/32) is worth reading.

If you have small enough datasets that the [R implemntation of DivNet](https://github.com/adw96/DivNet) can handle them, just use that instead...it will let you estimate the network structure of your data!

## License

Dual-licensed to be compatible with the Rust project.

Licensed under the Apache License, Version 2.0 or the MIT license, at your option. This program may not be copied, modified, or distributed except according to those terms.