use serde::Deserialize;
use std::path::PathBuf;

#[derive(Deserialize, Debug)]
pub struct Config {
    pub model: FitAitchisonConfig,
    pub io: IoConfig,
    pub misc: MiscConfig,
}

#[derive(Deserialize, Debug)]
pub struct FitAitchisonConfig {
    pub em_iter: usize,
    pub em_burn: usize,
    pub mc_iter: usize,
    pub mc_burn: usize,
    pub stepsize: f64,
    pub perturbation: f64,
    pub replicates: usize,
    pub base_taxa: usize,
}

#[derive(Deserialize, Debug)]
pub struct IoConfig {
    pub count_table: PathBuf,
    pub sample_data: PathBuf,
    pub output: PathBuf,
}

#[derive(Deserialize, Debug)]
pub struct MiscConfig {
    pub parallel_replicates: bool,
}
