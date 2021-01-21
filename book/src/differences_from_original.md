# Important differences from reference implementation

DivNet and `divnet-rs` are very similar and use the same algorithms.  However,
there are some differences.  I mention some of the important ones here.

## Input & Output files

The input and output files are different.  To see how to get files into
`divnet-rs` from R and then back in to R once `divnet-rs` is finished, see the
scripts in the [Lee example directory](./lee_data).

## Diagonal network only

The only network option is `diagonal`.  Here is an
[quote from Amy Willis](https://github.com/adw96/DivNet/issues/32#issuecomment-521727997)
that supports this decision:

> I would recommend network="diagonal" for a dataset of this size. This means
you're allowing overdispersion (compared to a plugin aka multinomial model) but
not a network structure. This isn't just about computational expense -- it's
about the reliability of the network estimates. Essentially estimating network
structure on 20k variables (taxa) with 50 samples with any kind of reliability
is going to be very challenging, and I don't think that it's worth doing here.
In our simulations we basically found that overdispersion contributes the bulk
of the variance to diversity estimation (i.e. overdispersion is more important
than network structure), so I don't think you are going to lose too much anyway.

The [whole issue](https://github.com/adw96/DivNet/issues/32) is worth reading.

If you have small enough datasets that the
[R implemntation of DivNet](https://github.com/adw96/DivNet) can handle them,
just use that instead...it will let you estimate the network structure of your 
data!

## Bootstrapping

Only the parametric bootstrap is available.  You cannot do the nonparametric
bootstrap.

## Monte-Carlo iterations and burn

The MC iterations must be even, and the MC burn must be 1/2 of the MC
iterations.  This allows me to do a little trick to cut the overall memory usage
by about 1/3.
