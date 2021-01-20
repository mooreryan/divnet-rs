extern crate blas;
extern crate divnet_rs;
extern crate lapack;
extern crate openblas_src;
#[macro_use]
extern crate log;

use divnet_rs::fit_aitchison;
use divnet_rs::fit_aitchison::parametric_bootstrap;
use divnet_rs::io;
use divnet_rs::opts::Opts;
use env_logger::Env;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::fs::File;
use std::io::{BufWriter, Write};
use structopt::StructOpt;

fn main() {
    // Set default log level to info.
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    let opts = Opts::from_args();

    debug!("opts: {:?}", &opts);

    let config = io::read_config(&opts.config);

    debug!("config: {:?}", &config);

    // TODO move these into the read_config function.
    assert!(config.model.mc_iter >= 2);
    assert!(config.model.mc_burn >= 1);
    assert_eq!(
        config.model.mc_iter - config.model.mc_burn,
        config.model.mc_burn,
        "MC burn must be 1/2 of MC iter"
    );

    // Todo switch to seedable rng
    let mut rng = ChaCha20Rng::seed_from_u64(config.misc.random_seed);

    let (sample_names, taxa_names, otu_table) = io::read_otu_table(&config.io.count_table);
    let (sample_names2, _variable_names, sample_data) =
        io::read_sample_data(&config.io.sample_data);

    if &sample_names[..] != &sample_names2[..] {
        panic!("sample name mismatch in input files!")
    }

    debug!("otu_table: {:?}", &otu_table.dim());
    debug!("sample_data: {:?}", &sample_data.dim());

    let outf = File::create(&config.io.output).expect("couldn't open output file for writing");
    let mut outf = BufWriter::new(outf);

    info!("Running fit_aitchison");
    let fa_result = fit_aitchison::fit_aitchison(&mut rng, &otu_table, &sample_data, &config.model);

    writeln!(&mut outf, "# this is replicate 0").unwrap();

    write!(&mut outf, "replicate,sample").unwrap();
    for taxa in taxa_names.iter() {
        write!(&mut outf, ",{}", &taxa).unwrap();
    }
    writeln!(&mut outf).unwrap();

    // the fa_result directly is the first replicate
    for sample_idx in 0..fa_result.fitted_z.ncols() {
        write!(&mut outf, "0,{}", sample_names[sample_idx]).unwrap();

        for taxa_idx in 0..fa_result.fitted_z.nrows() {
            write!(
                &mut outf,
                ",{}",
                fa_result.fitted_z.get(taxa_idx, sample_idx)
            )
            .unwrap();
        }

        writeln!(&mut outf).unwrap();
    }

    for bi in 0..config.model.replicates {
        info!("Doing replicate {}", bi + 1);

        let br = parametric_bootstrap(
            &mut rng,
            &fa_result.fitted_y,
            &fa_result.sigma_em,
            &otu_table,
            &sample_data,
            &config.model,
        );

        writeln!(&mut outf, "# this is replicate {}", bi + 1).unwrap();

        // the fa_result directly is the first replicate
        for sample_idx in 0..br.fitted_z.ncols() {
            write!(&mut outf, "{},{}", bi + 1, sample_names[sample_idx]).unwrap();

            for taxa_idx in 0..br.fitted_z.nrows() {
                write!(&mut outf, ",{}", br.fitted_z.get(taxa_idx, sample_idx)).unwrap();
            }

            writeln!(&mut outf).unwrap();
        }
    }

    info!("Done!  Output file is: {:?}", &config.io.output);
}
