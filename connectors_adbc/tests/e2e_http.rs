//! End-to-end HTTP/WebSocket tests against in-process server.
//! Run with: cargo test --test e2e_http -- --ignored
//! Starts the axum server on a random port and tests all endpoints.

use std::time::Duration;

use grok_connect_adbc::config::Config;
use grok_connect_adbc::providers::databricks::DatabricksProvider;
use grok_connect_adbc::providers::snowflake::SnowflakeProvider;
use grok_connect_adbc::providers::ProviderRegistry;
use grok_connect_adbc::state::{AppState, ConnectionPoolManager};

use axum::routing::{get, post};
use axum::Router;
use std::sync::Arc;
use tokio::net::TcpListener;

async fn start_server() -> String {
    let mut registry = ProviderRegistry::new();
    registry.register(Arc::new(DatabricksProvider::new()));
    registry.register(Arc::new(SnowflakeProvider::new()));

    let pool_manager = Arc::new(ConnectionPoolManager::new(5, 300));
    let config = Config::from_env();

    let state = AppState {
        registry: Arc::new(registry),
        pool_manager,
        config: Arc::new(config),
    };

    let app = Router::new()
        .route("/health", get(grok_connect_adbc::routes::health::health))
        .route("/info", get(grok_connect_adbc::routes::health::info))
        .route("/conn", get(grok_connect_adbc::routes::conn::list_providers))
        .route("/test", post(grok_connect_adbc::routes::test_connection::test_connection))
        .route("/query_socket", get(grok_connect_adbc::routes::query_socket::query_socket_upgrade))
        .route("/query", post(grok_connect_adbc::routes::query::query_stub))
        .route("/schema", post(grok_connect_adbc::routes::schema::schema))
        .route("/schemas", post(grok_connect_adbc::routes::schema::schemas))
        .route("/cancel", post(grok_connect_adbc::routes::cancel::cancel_stub))
        .with_state(state);

    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let url = format!("http://{}", addr);

    tokio::spawn(async move {
        axum::serve(listener, app).await.unwrap();
    });

    // Give server a moment to start
    tokio::time::sleep(Duration::from_millis(100)).await;
    url
}

#[tokio::test]
#[ignore]
async fn test_health_endpoint() {
    let url = start_server().await;
    let resp = reqwest::get(format!("{}/health", url)).await.unwrap();
    assert_eq!(resp.status(), 200);
    assert_eq!(resp.text().await.unwrap(), "OK");
}

#[tokio::test]
#[ignore]
async fn test_info_endpoint() {
    let url = start_server().await;
    let resp = reqwest::get(format!("{}/info", url)).await.unwrap();
    assert_eq!(resp.status(), 200);
    let json: serde_json::Value = resp.json().await.unwrap();
    assert_eq!(json["name"], "grok_connect_adbc");
    assert!(json["version"].is_string());
}

#[tokio::test]
#[ignore]
async fn test_conn_lists_providers() {
    let url = start_server().await;
    let resp = reqwest::get(format!("{}/conn", url)).await.unwrap();
    assert_eq!(resp.status(), 200);
    let descriptors: Vec<serde_json::Value> = resp.json().await.unwrap();

    let names: Vec<&str> = descriptors.iter()
        .filter_map(|d| d["type"].as_str())
        .collect();
    assert!(names.contains(&"DatabricksADBC"), "Should list DatabricksADBC");
    assert!(names.contains(&"SnowflakeADBC"), "Should list SnowflakeADBC");

    // Both should have canBrowseSchema=false
    for d in &descriptors {
        assert_eq!(d["canBrowseSchema"], false);
    }
}

#[tokio::test]
#[ignore]
async fn test_test_connection_bad_creds() {
    let url = start_server().await;
    let conn = serde_json::json!({
        "dataSource": "DatabricksADBC",
        "parameters": {
            "server": "invalid-host.databricks.com",
            "port": 443,
            "httpPath": "/sql/invalid"
        },
        "credentials": {
            "parameters": {
                "password": "invalid_token"
            }
        }
    });

    let client = reqwest::Client::new();
    let resp = client.post(format!("{}/test", url))
        .json(&conn)
        .send()
        .await
        .unwrap();
    // Should return 200 with error message, not 500
    assert_eq!(resp.status(), 200);
    let text = resp.text().await.unwrap();
    assert!(!text.contains("Connection available"), "Should fail with bad creds");
}

#[tokio::test]
#[ignore]
async fn test_query_socket_error() {
    let url = start_server().await;
    let ws_url = url.replace("http://", "ws://");

    use tokio_tungstenite::connect_async;
    use futures_util::{SinkExt, StreamExt};
    use tokio_tungstenite::tungstenite::Message;

    let (mut ws, _) = connect_async(format!("{}/query_socket", ws_url)).await.unwrap();
    let _ = ws.next().await; // CONNECTED

    // Send malformed QUERY
    ws.send(Message::Text("QUERY {invalid json}".to_string())).await.unwrap();

    let msg = ws.next().await.unwrap().unwrap();
    let text = msg.to_text().unwrap();
    assert!(text.starts_with("ERROR:"), "Expected ERROR message, got: {}", text);
}

#[tokio::test]
#[ignore]
async fn test_query_stub_returns_error() {
    let url = start_server().await;
    let client = reqwest::Client::new();
    let resp = client.post(format!("{}/query", url))
        .json(&serde_json::json!({}))
        .send()
        .await
        .unwrap();
    assert_eq!(resp.status(), 200);
    let result: serde_json::Value = resp.json().await.unwrap();
    assert_eq!(result["#type"], "DataQueryRunResult");
    assert!(result["errorMessage"].as_str().unwrap().contains("Not implemented"));
}

fn parse_framed(bytes: &[u8]) -> (serde_json::Value, Vec<u8>) {
    assert_eq!(&bytes[0..2], &0xA103u16.to_le_bytes(), "d42 uint8-list type code");
    let len = f64::from_le_bytes(bytes[2..10].try_into().unwrap()) as usize;
    let json = std::str::from_utf8(&bytes[10..10 + len]).unwrap();
    (serde_json::from_str(json).unwrap(), bytes[10 + len..].to_vec())
}

#[tokio::test]
#[ignore]
async fn test_schemas_unsupported_provider_returns_framed_error() {
    let url = start_server().await;
    let client = reqwest::Client::new();
    let resp = client.post(format!("{}/schemas", url))
        .json(&serde_json::json!({"dataSource": "SnowflakeADBC"}))
        .send()
        .await
        .unwrap();
    assert_eq!(resp.status(), 200);
    let bytes = resp.bytes().await.unwrap();
    let (result, payload) = parse_framed(&bytes);
    assert_eq!(result["#type"], "DataQueryRunResult");
    assert!(result["errorMessage"].as_str().unwrap().contains("Schema browsing not implemented"));
    assert!(payload.is_empty());
}

#[tokio::test]
#[ignore]
async fn test_schema_unknown_datasource_returns_framed_error() {
    let url = start_server().await;
    let client = reqwest::Client::new();
    let resp = client.post(format!("{}/schema", url))
        .json(&serde_json::json!({}))
        .send()
        .await
        .unwrap();
    assert_eq!(resp.status(), 200);
    let (result, _) = parse_framed(&resp.bytes().await.unwrap());
    assert!(result["errorMessage"].as_str().unwrap().contains("Unknown data source"));
}

#[tokio::test]
#[ignore]
async fn test_cancel_stub_returns_error() {
    let url = start_server().await;
    let client = reqwest::Client::new();
    let resp = client.post(format!("{}/cancel", url))
        .json(&serde_json::json!({}))
        .send()
        .await
        .unwrap();
    assert_eq!(resp.status(), 200);
    let result: serde_json::Value = resp.json().await.unwrap();
    assert!(result["errorMessage"].as_str().unwrap().contains("cancellation not implemented"));
}
