//! Integration tests for the ClickHouse provider (DLL loaded via the ADBC driver manager).
//! Run against a local ClickHouse (docker container `clickhouse`, localhost:8123,
//! default/clickhouse) with tables `bench_categorical` / `bench_longcat`:
//!   cargo test --test clickhouse_integration -- --ignored --test-threads=1
//! Host overridable via $CLICKHOUSE_HOST; driver path via $ADBC_CLICKHOUSE_DRIVER.

use arrow::datatypes::DataType;
use grok_connect_adbc::protocol::*;
use grok_connect_adbc::providers::clickhouse::ClickhouseProvider;
use grok_connect_adbc::providers::DbProvider;
use grok_connect_adbc::state::ConnectionPoolManager;
use std::collections::HashMap;

fn setup() -> (ClickhouseProvider, ConnectionPoolManager, DataConnection) {
    // Point the driver manager at the freshly-built DLL in lib/ (deterministic load).
    if std::env::var("ADBC_CLICKHOUSE_DRIVER").is_err() {
        let dll = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("lib")
            .join("adbc_clickhouse.dll");
        if dll.exists() {
            std::env::set_var("ADBC_CLICKHOUSE_DRIVER", dll);
        }
    }

    let host = std::env::var("CLICKHOUSE_HOST").unwrap_or_else(|_| "localhost".to_string());
    let port: i64 = std::env::var("CLICKHOUSE_PORT").ok().and_then(|s| s.parse().ok()).unwrap_or(8123);
    let password = std::env::var("CLICKHOUSE_PASSWORD").unwrap_or_else(|_| "clickhouse".to_string());
    let mut params = HashMap::new();
    params.insert("server".to_string(), serde_json::json!(host));
    params.insert("port".to_string(), serde_json::json!(port));

    let mut cred_params = HashMap::new();
    cred_params.insert("login".to_string(), serde_json::json!("default"));
    cred_params.insert("password".to_string(), serde_json::json!(password));

    let conn = DataConnection {
        id: "test".to_string(),
        data_source: "ClickHouseADBC".to_string(),
        connection_string: None,
        credentials: Some(Credentials { parameters: cred_params }),
        parameters: params,
        tags: HashMap::new(),
    };

    (ClickhouseProvider::new(), ConnectionPoolManager::new(5, 300), conn)
}

fn make_call(conn: DataConnection, query: &str) -> FuncCall {
    FuncCall {
        id: "test".to_string(),
        func: DataQuery {
            type_tag: "DataQuery".to_string(),
            id: String::new(),
            name: String::new(),
            query: query.to_string(),
            connection_id: String::new(),
            connection: conn,
            params: vec![],
            options: HashMap::new(),
            aux: HashMap::new(),
        },
        options: HashMap::new(),
        parameter_values: HashMap::new(),
        aux: HashMap::new(),
        log: None,
    }
}

#[test]
#[ignore]
fn test_clickhouse_connectivity() {
    let (provider, pool_mgr, conn) = setup();
    provider.test_connection(&conn, &pool_mgr).expect("Connection test failed");
}

#[test]
#[ignore]
fn test_clickhouse_select_literal() {
    let (provider, pool_mgr, conn) = setup();
    let call = make_call(conn, "SELECT 42 AS answer");
    let batches = provider.execute_query(&call, &pool_mgr).expect("Query failed");
    assert_eq!(batches.iter().map(|b| b.num_rows()).sum::<usize>(), 1);
}

/// The whole point: LowCardinality flows through the connector as an Arrow Dictionary
/// (the provider issues `SET output_format_arrow_low_cardinality_as_dictionary=1` on connect),
/// while a plain String column stays Utf8 — all read via arrow 57 across the C ABI.
#[test]
#[ignore]
fn test_clickhouse_lowcardinality_is_dictionary() {
    let (provider, pool_mgr, conn) = setup();
    let call = make_call(conn, "SELECT status, uuid_col FROM bench_categorical LIMIT 1000");
    let batches = provider.execute_query(&call, &pool_mgr).expect("Query failed");
    let total: usize = batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(total, 1000);

    let schema = batches[0].schema();
    let status_type = schema.field(schema.index_of("status").unwrap()).data_type().clone();
    let uuid_type = schema.field(schema.index_of("uuid_col").unwrap()).data_type().clone();

    assert!(
        matches!(status_type, DataType::Dictionary(_, _)),
        "LowCardinality(String) `status` should arrive as Dictionary, got {status_type:?}"
    );
    assert_eq!(
        uuid_type,
        DataType::Utf8,
        "plain String `uuid_col` should stay Utf8, got {uuid_type:?}"
    );
}
