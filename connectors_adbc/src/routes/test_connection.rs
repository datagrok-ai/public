use axum::extract::State;
use axum::response::IntoResponse;
use axum::http::StatusCode;
use axum::Json;
use crate::state::AppState;
use crate::protocol::DataConnection;
use tracing::{info, error};

/// POST /test — test database connection.
/// Expects DataConnection JSON body. Returns "Connection available" or error message.
pub async fn test_connection(
    State(state): State<AppState>,
    Json(conn): Json<DataConnection>,
) -> impl IntoResponse {
    let data_source = conn.data_source.clone();
    info!(provider = %data_source, "Testing connection");

    let provider = match state.registry.get(&data_source) {
        Some(p) => p.clone(),
        None => {
            let msg = format!("Unknown data source: {}", data_source);
            error!(%msg);
            return (StatusCode::BAD_REQUEST, msg);
        }
    };

    let pool_mgr = state.pool_manager.clone();

    // Run sync ADBC call on blocking thread
    let result = tokio::task::spawn_blocking(move || {
        provider.test_connection(&conn, &pool_mgr)
    }).await;

    match result {
        // "ok" is the success token Datlas/the client expect (DataProvider.CONN_AVAILABLE),
        // matching Java grok_connect. Returning any other string surfaces as a failed test.
        Ok(Ok(())) => (StatusCode::OK, "ok".to_string()),
        Ok(Err(e)) => {
            let msg = format!("{:#}", e);
            error!(provider = %data_source, error = %msg, "Connection test failed");
            (StatusCode::OK, msg)
        }
        Err(e) => {
            let msg = format!("Internal error: {}", e);
            error!(%msg);
            (StatusCode::INTERNAL_SERVER_ERROR, msg)
        }
    }
}
