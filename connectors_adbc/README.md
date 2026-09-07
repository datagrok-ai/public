# Grok Connect ADBC

A Rust implementation of the Grok Connect protocol built on [ADBC](https://arrow.apache.org/adbc/)
(Arrow Database Connectivity). It serves the same REST/WebSocket endpoints as the Java
[grok_connect](../connectors/README.md), so the platform routes queries to it transparently —
by default it is registered as an additional endpoint alongside the Java one, and each
provider type is served by exactly one endpoint.

Data flows through the service as Arrow record batches end-to-end: the ADBC driver returns
Arrow, and the service streams it to the platform without row-by-row conversion.

## Supported providers

| Provider | Driver |
|----------|--------|
| Snowflake | `adbc_driver_snowflake` dynamic library, loaded at runtime |
| Databricks | Native Rust ADBC driver |
| BigQuery | C-ABI dynamic library, loaded via the ADBC driver manager |
| ClickHouse | C-ABI dynamic library, loaded via the ADBC driver manager |

Dynamic drivers are fetched into `lib/` with `scripts/fetch_adbc_driver.sh <provider>`
(`fetch_adbc_driver.ps1` on Windows), or their location can be overridden with
`ADBC_SNOWFLAKE_DRIVER` / `ADBC_BIGQUERY_DRIVER` / `ADBC_CLICKHOUSE_DRIVER`.

## Build and run

```bash
cargo build --release
./target/release/grok_connect_adbc
```

Or with Docker:

```bash
docker build -t datagrok/grok_connect_adbc .
docker run -p 1234:1234 datagrok/grok_connect_adbc
```

Configuration is via environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `GROK_CONNECT_PORT` | `1234` | HTTP/WebSocket listen port |
| `POOL_MAX_SIZE` | `30` | Max pooled connections per target |
| `POOL_IDLE_TIMEOUT_SECS` | `300` | Idle connection eviction timeout |
| `POOL_EVICTION_INTERVAL_SECS` | `60` | Pool eviction sweep interval |
| `RUST_LOG` | | Log filter, e.g. `grok_connect_adbc=info` |

Database credentials are never configured on the service — they arrive per request
in the `FuncCall` payload, same as with the Java grok_connect.

## Tests

```bash
cargo test                    # unit + protocol tests
./scripts/run_all_tests.sh    # full suite, including Docker-based integration tests
```

Integration tests that need live database credentials are skipped unless the
corresponding environment variables (e.g. `SNOWFLAKE_CREDS_PATH`) are set — see
`scripts/e2e_smoke.sh`.
