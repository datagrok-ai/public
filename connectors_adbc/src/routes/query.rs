use axum::response::Json;
use crate::protocol::DataQueryRunResult;

/// POST /query — stub endpoint. Returns structured error directing to WebSocket endpoint.
/// Keeps API surface consistent with Java grok_connect so Datlas doesn't get 404s.
pub async fn query_stub() -> Json<DataQueryRunResult> {
    Json(DataQueryRunResult::error(
        "Not implemented in grok_connect_adbc. Use WebSocket /query_socket endpoint for streaming Arrow IPC results."
    ))
}
