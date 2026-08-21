use axum::response::Json;
use crate::protocol::DataQueryRunResult;

/// POST /cancel — stub. Query cancellation not implemented in PoC.
pub async fn cancel_stub() -> Json<DataQueryRunResult> {
    Json(DataQueryRunResult::error(
        "Query cancellation not implemented in grok_connect_adbc PoC."
    ))
}
