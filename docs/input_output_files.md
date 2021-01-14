# Input & Output Files

There are scripts to help you generate the input files for `divnet-rs` from R. 
I recommend using them rather than making input files by hand!  There are also
scripts to help you re-import the output of `divnet-rs` back into R.

- [ ] TODO upload the helper scripts.

## Input files

The `divnet-rs` program takes two input files:  a count table and a file with
sample data.  Both should be CSV files.

### Count table

The count table (aka OTU table, aka ASV table) is a taxa-by-sample
representation of your experiment.

Here is an example:

```csv
taxa,s1,s2,s3
t1,200,1,2
t2,210,2,1
t3,180,2,1
t4,1,230,235
t5,2,220,215
```

This data set has five taxa (`t1`, `t2`, `t3`, `t4`, `t5`) and three samples
(`s1`, `s2`, `s3`).  Note that these can be named whatever you want.

The `taxa` specifier is ignored and so you can write whatever you want there.
E.g., if you have amplicon sequence variants, you could put `asv` there instead
of `taxa`.

The values are counts, so they should be positive integers only.

- [ ] TODO check if spaces in names are allowed.
- [ ] TODO check if names have to be unique.
- [ ] TODO check if it works with relative abundances.

### Sample data

The sample data file is a little weird looking, I will admit.  It is basically
the output of the R function `model.matrix`.  It converts your dummy variables
to 0 and 1.

#### A simple example 

In this case, there is only one covariate: `snazzy`.  Here, I have labeled it as
`snazzyyes` indicating that samples with a `1` are snazzy (i.e., positive for
condition `snazzy`) and samples with a `0` are not snazzy (negative for
condition `snazzy`).  So `s1` is snazzy, but `s2` and `s3` are NOT snazzy.

```csv
sample,snazzyyes
s1,1
s2,0
s3,0
```

- [ ] TODO can you include continuous variables or just categorical data?

#### Another example

Here are a couple of lines from the sample data file for the Lee dataset that's
included in DivNet.

```csv
sample,charbiofilm,charcarbonate,charglassy,charwater
BW1,0,0,0,1
BW2,0,0,0,1
R10,0,0,1,0
R11,0,0,1,0
```

As you can see, the variable of interest is `char`.  It has the following
columns:

- `charbiofilm` (`1` for yes, it's a biofilm sample, `0` for no it is not)
- `charcarbonate` (`1` for yes, it's a carbonate sample, `0` for no it is not)
- `charglassy` (`1` for yes, it's a glassy sample, `0` for no it is not)
- `charwater` (`1` for yes, it's a water sample, `0` for no it is not)

Now the Lee data has a fifth category, `alered`.  It is not listed here as that
is the way the `model.matrix` dummy encoding works.  You don't need a column for
it, (and if you do include it in your dummy encoding things can get wonky) any
sample with a `0` in all the colunms is an `altered` sample.

## Output files

Here is the output file you get if you run the example files in `<source
root>/test_files/small`.  (The little ones you see above!)

```csv
# this is replicate 0
replicate,sample,t1,t2,t3,t4,t5
0,s1,0.33855749713477606,0.34823818972504217,0.30865008153567924,0.0018295805613595625,0.0027246510431431564
0,s2,0.0030379849794267364,0.0028353087370401467,0.0033011806108468353,0.490330204063212,0.5004953216094743
0,s3,0.0030379849794267364,0.0028353087370401467,0.0033011806108468353,0.490330204063212,0.5004953216094743
# this is replicate 1
1,s1,0.3627204546824228,0.36468668498369233,0.26626637816414933,0.00008013857010885896,0.006246343599626745
1,s2,0.0006637323562423841,0.0027111558898110805,0.003702116838977756,0.5432176368099867,0.44970535810498213
1,s3,0.0006637323562423841,0.0027111558898110805,0.003702116838977756,0.5432176368099867,0.44970535810498213
# this is replicate 2
2,s1,0.5294663181856507,0.2428125790405528,0.22595059438495904,0.0016861984655594231,0.00008430992327797039
2,s2,0.003824488773323916,0.003824488773323916,0.005408643892378329,0.5034205045920187,0.48352187396895513
2,s3,0.003824488773323916,0.003824488773323916,0.005408643892378329,0.5034205045920187,0.48352187396895513
# this is replicate 3
3,s1,0.5901081082784421,0.22421302607955013,0.17869297960014005,0.0015961885594974962,0.005389697482370069
3,s2,0.0030619926513927457,0.0006994389587677033,0.004747225069582607,0.5271925248691112,0.4642988184511458
3,s3,0.0030619926513927457,0.0006994389587677033,0.004747225069582607,0.5271925248691112,0.4642988184511458
# this is replicate 4
4,s1,0.42087439504861085,0.3408217058601948,0.22533144275449316,0.0035274481524021776,0.009445008184299079
4,s2,0.004061119432165891,0.0007170056373956505,0.0007467384699807625,0.5003428103976993,0.4941323260627584
4,s3,0.004061119432165891,0.0007170056373956505,0.0007467384699807625,0.5003428103976993,0.4941323260627584
# this is replicate 5
5,s1,0.20558302906768744,0.49136400424687354,0.29390176182754923,0.005463507248674078,0.0036876976092156933
5,s2,0.000702351053334611,0.0024864420846098466,0.0032837333411371906,0.43124409424160054,0.5622833792793178
5,s3,0.000702351053334611,0.0024864420846098466,0.0032837333411371906,0.43124409424160054,0.5622833792793178
```

Again, not all that nice for human consumption, but it will be nice and easy to
parse in R.  Check out the helper scripts for that!

- [ ] TODO upload helper scripts.
