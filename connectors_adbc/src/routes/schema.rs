use axum::extract::{Query, State};
use axum::http::header;
use axum::response::{IntoResponse, Response};
use axum::Json;
use crate::protocol::{frame_with_header, DataConnection, DataQueryRunResult, FuncCall};
use crate::state::AppState;
use std::collections::HashMap;
use tracing::{error, info};

pub async fn schemas(
    State(state): State<AppState>,
    Json(conn): Json<DataConnection>,
) -> Response {
    let sql = state.registry.get(&conn.data_source).and_then(|p| p.schemas_sql());
    run_schema_query(state, conn, sql).await
}

pub async fn schema(
    State(state): State<AppState>,
    Query(q): Query<HashMap<String, String>>,
    Json(conn): Json<DataConnection>,
) -> Response {
    let include_key_info = q.get("includeKeyInfo").map(String::as_str) == Some("true");
    let schema = conn.get_str("schema").map(str::to_string);
    let table = conn.get_str("table").map(str::to_string);
    let sql = state
        .registry
        .get(&conn.data_source)
        .and_then(|p| p.schema_sql(schema.as_deref(), table.as_deref(), include_key_info));
    run_schema_query(state, conn, sql).await
}

async fn run_schema_query(state: AppState, conn: DataConnection, sql: Option<String>) -> Response {
    let data_source = conn.data_source.clone();
    let provider = match state.registry.get(&data_source) {
        Some(p) => p.clone(),
        None => return framed_error(format!("Unknown data source: {}", data_source)),
    };
    let Some(sql) = sql else {
        return framed_error(format!(
            "Schema browsing not implemented for {} (canBrowseSchema=false).",
            data_source
        ));
    };

    info!(provider = %data_source, sql = %sql, "Schema request");
    let pool_mgr = state.pool_manager.clone();
    let call = FuncCall::for_raw_sql(sql, conn);
    let started = std::time::Instant::now();
    let result =
        tokio::task::spawn_blocking(move || provider.execute_query(&call, &pool_mgr)).await;

    let batches = match result {
        Ok(Ok(batches)) => batches,
        Ok(Err(e)) => return framed_error(format!("{:#}", e)),
        Err(e) => return framed_error(format!("Internal error: {}", e)),
    };
    let rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    let columns = batches.first().map(|b| b.num_columns()).unwrap_or(0);
    match crate::arrow_ipc::record_batches_to_ipc(&batches) {
        Ok(ipc) => {
            let mut run = DataQueryRunResult::success(
                columns as i32, rows as i32, started.elapsed().as_secs_f64());
            run.blob_length = ipc.len() as i64;
            framed(&run, &ipc)
        }
        Err(e) => framed_error(format!("{:#}", e)),
    }
}

fn framed(run: &DataQueryRunResult, payload: &[u8]) -> Response {
    let json = serde_json::to_string(run).unwrap();
    (
        [(header::CONTENT_TYPE, "application/octet-stream")],
        frame_with_header(&json, payload),
    )
        .into_response()
}

fn framed_error(message: String) -> Response {
    error!(%message, "Schema request failed");
    framed(&DataQueryRunResult::error(message), &[])
}
