extern crate blas;
extern crate divnet_rs;
extern crate lapack;
extern crate openblas_src;
#[macro_use]
extern crate log;

use bincode::serialize_into;
use divnet_rs::fit_aitchison;
use divnet_rs::io;
use divnet_rs::opts::Opts;
use env_logger::Env;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use structopt::StructOpt;

fn main() {
    // Set default log level to info.
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    let opts = Opts::from_args();

    debug!("opts: {:?}", &opts);

    let mut config = io::read_config(&opts.config);

    debug!("config: {:?}", &config);

    let fa_result_dat = env::var("DN_MODEL_DAT").unwrap_or("/tmp/dn_model".to_string());
    let fa_result_csv = env::var("DN_MODEL_CSV").unwrap_or("/tmp/dn_model".to_string());

    config.io.output = PathBuf::from(fa_result_csv);

    // TODO move these into the read_config function.
    assert!(config.model.mc_iter >= 2);
    assert!(config.model.mc_burn >= 1);
    assert_eq!(
        config.model.mc_iter - config.model.mc_burn,
        config.model.mc_burn,
        "MC burn must be 1/2 of MC iter"
    );

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

    let mut f = BufWriter::new(File::create(fa_result_dat).unwrap());
    serialize_into(&mut f, &fa_result).unwrap();

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

    info!("Done!  Output file is: {:?}", &config.io.output);
}
