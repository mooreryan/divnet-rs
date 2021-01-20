# Logging

`divnet-rs` logs various things to standard error during a run.

## Log levels

There are multiple log levels.  Here they are from least important to most important.

- `trace`
- `debug`
- `info`
- `warn`
- `error`

By default `divnet-rs` only prints `info`, `warn`, and `error` messages.

If you want, you can turn on the lower level messages using the `RUST_LOG` [environmental variable](https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-linux) like this:

```
RUST_LOG=debug divnet-rs config.toml
```

That would tell `divnet-rs` to print all `debug` messages plus any messages that are of higher importance (`info`, `warn`, and `error`).

If you want fewer messages, you could use `RUST_LOG=warn` or `RUST_LOG=error`.

## Recommendation

Generally the default log level is fine.  It's also nice since you don't have to set any environmental variables!

If you have very large input files and want to see the most logging, for example to try and estimate the time remaining, you could turn on mode messages with `RUST_LOG=trace` or `RUST_LOG=trace`.  It will put a lot of output though.

## Setting config level and number of threads

You set the number of threads for OpenBLAS with an environment variable.  If you want to set the log level and the OpenBLAS threads, you can do it like this:

```shell
OPENBLAS_NUM_THREADS=1 RUST_LOG=debug divnet-rs /path/to/config.toml
```

See [Config Files](./config_files.md) for more info on the OpenBLAS thread options.
