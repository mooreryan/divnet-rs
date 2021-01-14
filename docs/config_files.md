# Configuration Files

You control `divnet-rs` using a [TOML](https://toml.io/) configuration file
(config file).

Here is an example of one with lots of comments to help you see what's going on.
 Feel free to copy this file and use it as a model for your own datasets!

```toml
[model]
# Number of expectation maximization (EM) iterations
em_iter = 6

# Number of EM iterations to burn (i.e., throw out)
em_burn = 3

# Number of Monte-Carlo (MC) iterations.
mc_iter = 500

# Number of MC iterations to burn
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
# If true, then do the replicates in parallel with rayon.  As of now, you should
# probably just leave this as false.  On my four core test macine, setting it to
# true saves a BIT of walltime at the expense of double the CPU time.
parallel_replicates = false
```

As you can see it is broken up in to three sections (model, io, and misc) each
controlling a different aspect of the software.  The `[model]` section contains
config options for the actual model DivNet uses to estimate diversity.  The
`[io]` section deals with specifying input and output files.  Finally, the
`[misc]` section contains any random options that don't fit in any other
category.  Currently the only option in this category is `parallel_replicates`. 
Note that while this sounds cool, it actually gives worse performance right now,
so you should just always set it to false as you see above.

## Defaults

Currently there *are no default values*.  This means that every config file you
write has to explicitly list all of the options you see in the file above.  I
would like to eventually change this so you only have to specify values that
differ from the defaults, but as of now, you will have to be explicit :)

That said, DivNet does have some defaults for the tuning of the algorithm.  Here
is how you would set up the `[model]` section to use the defaults.

One interesting thing is that `divnet-rs` is *MUCH* faster than the R
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