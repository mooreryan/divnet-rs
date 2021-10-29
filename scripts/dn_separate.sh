DN_MODEL_DAT=dn_model DN_MODEL_CSV=dn_model.csv ./target/release/fit_model ./test_files/small/config_small.toml

# Each of these could be run on separate nodes, or as separate scheduled jobs (eg torque, slurm, etc).
#
# You must change the replicate, seed, and outfiles by hand!
DN_MODEL_DAT=dn_model DN_REPLICATE=1 DN_SEED=1 DN_BOOT_CSV=dn_boot_1.csv  ./target/release/bootstrap ./test_files/small/config_small.toml
DN_MODEL_DAT=dn_model DN_REPLICATE=2 DN_SEED=2 DN_BOOT_CSV=dn_boot_2.csv  ./target/release/bootstrap ./test_files/small/config_small.toml
DN_MODEL_DAT=dn_model DN_REPLICATE=3 DN_SEED=3 DN_BOOT_CSV=dn_boot_3.csv  ./target/release/bootstrap ./test_files/small/config_small.toml

# Concatenate the outfiles
cat dn_model.csv dn_boot_*.csv > dn_out.csv

# Clean up
rm dn_model dn_model.csv dn_boot_*.csv

