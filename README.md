# divnet-rs

It's [DivNet](https://github.com/adw96/DivNet) with fewer features, but also *way* faster and more memory efficient!!

Unfortunately, it is also more complicated to use at the moment, so if your data sets are smaller, I recommend using the [original implementation](https://github.com/adw96/DivNet).

Please see the [GitHub repository for the reference implementation](https://github.com/adw96/DivNet).  It has lots of good info about the theory behind DivNet!

## Installation & Usage

For installation and usage instructions, please see the [divnet-rs book](https://mooreryan.github.io/divnet-rs-book/).

### Default branch is now main

The default branch is now `main`, and the `master` branch no longer exists.

If you have a local clone using `master` as the default branch, you can update it by running the following commands.

```
git branch -m master main
git fetch origin
git branch -u origin/main main
git remote set-head origin -a
```

Optionally, run the following command to remove tracking references to the old branch name.

```
git remote prune origin
```

## Citation

If you use DivNet, please support the authors and [cite their manuscript](https://doi.org/10.1093/biostatistics/kxaa015)!

## License

Licensed under the Apache License, Version 2.0 or the MIT license, at your option. This program may not be copied, modified, or distributed except according to those terms.
