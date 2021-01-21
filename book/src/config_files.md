# Configuration Files

You control `divnet-rs` using a [TOML](https://toml.io/) configuration file
(config file).

Here is an example of one with lots of comments to help you see what's going on.
 Feel free to copy this file and use it as a model for your own datasets!

```toml
[model]
# Number of expectation maximization (EM) iterations
em_iter = 6

# Number of EM iterations to burn (i.e., throw out). Unlike mc_burn, this does
# NOT have to be em_iter / 2.
em_burn = 3

# Number of Monte-Carlo (MC) iterations.  Must be even (see mc_burn for
# details).
mc_iter = 500

# Number of MC iterations to burn. It must be mc_iter / 2.  If not, the program
# will abort.
mc_burn = 250

# Variance used for MH samples
stepsize = 0.01

# Perterbation magnitude for zero values when calculating logratios. (Any zero
# value will be replaced with this.) 
perturbation = 0.05

# Number of bootstrap iterations for estimating the variance.
replicates = 5

# The "base" taxa to use in the logratios.  The number represents the (0-based)
# index of the taxon you want.  So 0 here means the first taxon, 1 means the
# second, and so on.  Ideally it is a taxa observed in all samples.  That's not
# likely though, so try a couple of highly abundant taxa to confirm the results
# are robust to the taxon choice. 
base_taxa = 0

[io]
# If you use relative paths here, they will be interpreted relative to the
# directory from which you run the divnet-rs command.  If things get weird with
# that, or you're using a job scheduler like slurm or torque, just specify the
# absolute path here instead.

# The path to your count/asv/otu table
count_table = "./test_files/small/small_count_table.csv"

# The path to your sample data
sample_data = "./test_files/small/small_sample_data.csv"

# The path to the file that divnet-rs will generate
output = "./test_files/small/small_divnet_output.csv"

[misc]
# An unsigned 64 bit integer used to seed the random number generator.
random_seed = 0
```

As you can see it is broken up in to three sections (model, io, and misc) each
controlling a different aspect of the software.  The `[model]` section contains
config options for the actual model DivNet uses to estimate diversity.  The
`[io]` section deals with specifying input and output files.  Finally, the
`[misc]` section contains any options that don't fit in any other
category.

**IMPORTANT**:  The program will abort unless `mc_burn == mc_iter / 2`.  This is
*allows me to do some trickery to reduce the overall memory usage by ~1/3.

### Number of threads

One thing that isn't in the config files, but something you will want to do is to run `divnet-rs` with the `OPENBLAS_NUM_THREADS=1` environmental variable set.  This will ensure that OpenBLAS is only using a single thread.  In all my tests, `divnet-rs` will be between 25-50% faster with 1 thread for OpenBLAS.  You run it like this:

```
OPENBLAS_NUM_THREADS=1 divnet-rs /path/to/config.toml
```

Or you could just set that in your system config files.

See [Logging](./logging.md) for how to combine this variable with the logging environmental variables.

## Defaults

Currently there *are no default values*.  This means that every config file you
write has to explicitly list all of the options you see in the file above.  I
would like to eventually change this so you only have to specify values that
differ from the defaults, but as of now, you will have to be explicit :)

That said, DivNet does have some defaults for the tuning of the algorithm.  Here
is how you would set up the `[model]` section to use the defaults.

Note that `divnet-rs` is ~20x faster than the R
implementation.  This means that depending on the size of your dataset, you
could crank `em_*`, `mc_*`, and `replicates` options up really high and see what
happens.  I haven't really tried this out yet, but it might be interesting!

### DivNet 'fast' default

```toml
[model]
em_iter = 6
em_burn = 3
mc_iter = 500
mc_burn = 250
stepsize = 0.01
perturbation = 0.05
replicates = 5
base_taxa = 0
```

### DivNet 'careful' default

```toml
[model]
em_iter = 10
em_burn = 5
mc_iter = 1000
mc_burn = 500
stepsize = 0.01
perturbation = 0.05
replicates = 5
base_taxa = 0
```
