use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(name = "divnet-rs", about = "A fast implementation of DivNet")]
pub struct Opts {
    /// Path to config-lee-phylum.toml
    #[structopt(parse(from_os_str))]
    pub config: PathBuf,
}
