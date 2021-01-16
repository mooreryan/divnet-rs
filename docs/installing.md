# Installing divnet-rs

I have successfully installed and run `divnet-rs` on both Debian Linux and
MacOSX.  I haven't tested this out on Windows at all :)

## Dependencies

`divnet-rs` is written in [Rust](https://www.rust-lang.org/).  If you do not
have Rust, you must [install it first](https://www.rust-lang.org/tools/install).

`divnet-rs` depends on BLAS and LAPACK for fast numerical computation.  While there are lots of different implementations, `divnet-rs` uses [OpenBLAS](https://www.openblas.net/).  To install OpenBLAS, you will need a C compiler and a Fortran compiler.  I recommend [GCC, the GNU Compiler Collection](https://gcc.gnu.org/) for this.  With the GCC toolchain, I was able to install on both Linux and Mac with no issues.

*Note:  If you already have gcc, and you're using a Mac, you should check that you have an actual gcc program.  Often, you will use `gcc` from the command line, but it will actually be [clang](https://clang.llvm.org/).  You can check by running `gcc --vertion` and inspecting the output.  I was **not** able to compile `divnet-rs` using clang.*

### Linux

If you don't already have `gcc` you could install it from your distro's package manager e.g., `sudo apt-get install build-essentials` or something similar.  This will get you both `gcc` and `gfortran`. Also, you can check out the [official gcc install docs](https://gcc.gnu.org/install/).  

Additionally you will need to install OpenBLAS.  You can install it [from source](TOOD), but many distros also have it in their package managers. 

You may need to adjust the `LIBRARY_PATH` and `LD_LIBRARY_PATH` environmental variables depending on how you installed the above software.  See below for more info about that. 

### Mac

I used [Homebrew](https://brew.sh/) to install GCC and OpenBLAS like so:

```
brew install gcc
brew install openblas
```

Given the way Homebrew installs these packages, I needed to tweak some environment variables so that Cargo and the Rust compiler could actually link to the libraries.

When I installed with Homebrew, the OpenBLAS libraries were in `/usr/local/opt/openblas/lib` and the GCC libraries were in `/usr/local/lib/gcc/10`.  So I needed to set the `LIBRARY_PATH` and `LD_LIBRARY_PATH` environmental variables like so:

```shell
# For building
LIBRARY_PATH="/usr/local/opt/openblas/lib"
LIBRARY_PATH="${LIBRARY_PATH}:/usr/local/lib/gcc/10"
export LIBRARY_PATH

# For running
LD_LIBRARY_PATH="/usr/local/opt/openblas/lib"
LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib/gcc/10"
export LD_LIBRARY_PATH
```

Note that if you already have them set, you will need to append rater than overwrite them (as I did when setting the OpenBLAS paths).

You should put that in either `~/.profile`, `~/.bash_profile`, `~/.zprofile`, or whatever place you keep your shell startup configuration.

## Installing from source

There are two options.  You can either download a tagged release (recommended) or follow the master branch if you want the latest updates.

If you want to track the master branch, clone the git repository like so:

```
git clone https://github.com/mooreryan/divnet-rs.git
```

However, I recommend that you download a
[release](https://github.com/mooreryan/divnet-rs/releases) as they are versioned.

Once you have the source code, `cd` in to the root of the source directory and
use `cargo` to compile the program. (Cargo will be installed when you install Rust.)  Here is an example:

```
cd divnet-rs
cargo build --release
```

Note that this can take a while, but you will *not* have to compile the program
each time you use it!

You should now have a compiled binary in `./target/release/divnet-rs`.  You will
probably want move or symlink the executable binary
(`./target/release/divnet-rs`) somewhere on your path.  For example, if you want
to create a symlink to `~/bin`, you could run something like this:

```
ln -s $(pwd)/target/release/divnet-rs $HOME/bin/divnet-rs
```

You can check that you have access to the `divnet-rs` binary with the `which`
command:

```
which divnet-rs #=> /home/ryan/bin/divnet-rs
```

If you see a path, then the linking worked!  If not, then you can still just use
the program from its current directory (`target/release/divnet-rs`) or try and
[figure out how to get it on your path](https://gist.github.com/nex3/c395b2f8fd4b02068be37c961301caa7).

## Testing your installation

The `divnet-rs` repository includes some handy test files to check that
everything is working correctly.

Try it out!

```
divnet-rs ./test_files/small/config_small.toml
```

If all goes well, you should see some logging output that ends with a line like
this:

```
[2021-01-14T02:59:49Z INFO  divnet_rs] Done!  Output file is: "./test_files/small/small_divnet_output.csv"
```

Additionally there are some automated tests that you can run with `cargo test`.