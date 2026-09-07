use axum::extract::State;
use axum::extract::ws::{Message, WebSocket, WebSocketUpgrade};
use axum::response::IntoResponse;
use crate::dict_encode::{maybe_encode_batch, DictEncodeConfig};
use crate::protocol::{DataQueryRunResult, FuncCall};
use crate::state::AppState;
use tracing::{info, error, warn};

/// Buffer of 32 lets the ADBC producer race ahead of the WebSocket sender
/// (which is rate-limited by Datlas round-trips). Otherwise the slow consumer
/// would back-pressure CloudFetch downloads.
const BATCH_CHANNEL_SIZE: usize = 32;

/// Buffer of 2 gives perfect pipelining (drain produces chunk N+1 while main
/// sends chunk N) with one slot of jitter tolerance, without bloating RAM.
const CHUNK_CHANNEL_SIZE: usize = 2;

/// Chunking strategy — mirrors Java `grok_connect`'s `QueryManager`:
///
/// - Constants match Java: 10 MB adaptive byte target, [100, 100_000] row
///   bounds, small initial-chunk probe for fast time-to-first-row.
/// - On localhost both hops are free so bigger chunks win, but in production
///   each chunk hops grok_connect_adbc → Datlas → browser, and oversized chunks
///   head-of-line-block the pipeline. A 10 MB cap strikes the Java-compatible
///   balance: big enough to amortise per-chunk WS ceremony overhead, small
///   enough that Datlas can relay chunk N while we produce chunk N+1.
/// - The initial chunk is deliberately tiny so the client receives *something*
///   almost immediately (Java uses `INIT_FETCH_SIZE = 100` rows; we use a
///   ~200 KB byte probe since we accumulate driver-supplied Arrow batches
///   instead of requesting row counts from JDBC).
/// - Both targets are overridable per-query via `connectFetchSize` /
///   `initConnectFetchSize` in the FuncCall options (matching Java), and
///   per-process via `CHUNK_TARGET_BYTES` / `INIT_CHUNK_TARGET_BYTES` env vars.
const MAX_CHUNK_SIZE_BYTES: usize = 10_000_000;
const INIT_CHUNK_TARGET_BYTES_DEFAULT: usize = 200_000;
const MAX_FETCH_ROWS: usize = 100_000;
const MIN_FETCH_ROWS: usize = 100;

/// Chunking strategy for a single query. Resolved from env vars and per-query
/// FuncCall options before the drain task starts.
#[derive(Debug, Clone, Copy)]
struct ChunkStrategy {
    /// Target bytes for the first chunk (fast TTFR probe).
    init_target_bytes: usize,
    /// Target bytes for all subsequent chunks.
    target_bytes: usize,
    /// Hard cap on rows per chunk regardless of byte size.
    max_rows: usize,
    /// Optional exact row count per chunk (from FuncCall option
    /// `connectFetchSize = <int>`); overrides `target_bytes` when set.
    exact_rows: Option<usize>,
    /// Optional exact row count for the first chunk (from
    /// `initConnectFetchSize = <int>`); overrides `init_target_bytes`.
    init_exact_rows: Option<usize>,
    /// Dictionary-encode low-cardinality `Utf8` columns before IPC serialisation.
    /// `None` disables the pass entirely. Default = `Some(DictEncodeConfig::default())`
    /// unless overridden by `DICT_ENCODE_STRINGS` env or `meta.dictEncodeStrings`
    /// per-query option.
    dict_encode: Option<DictEncodeConfig>,
}

impl ChunkStrategy {
    fn from_env_and_call(call: &FuncCall) -> Self {
        let env_target = std::env::var("CHUNK_TARGET_BYTES")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(MAX_CHUNK_SIZE_BYTES);
        let env_init = std::env::var("INIT_CHUNK_TARGET_BYTES")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(INIT_CHUNK_TARGET_BYTES_DEFAULT);

        let dict_encode = std::env::var("DICT_ENCODE_STRINGS")
            .ok()
            .map(|v| DictEncodeConfig::from_str_option(&v))
            .unwrap_or(Some(DictEncodeConfig::default()));

        let mut strat = ChunkStrategy {
            init_target_bytes: env_init,
            target_bytes: env_target,
            max_rows: MAX_FETCH_ROWS,
            exact_rows: None,
            init_exact_rows: None,
            dict_encode,
        };

        // Per-query overrides from FuncCall.options, same keys as Java:
        //   connectFetchSize:      "N MB" → byte target; "<int>" → row count.
        //   initConnectFetchSize:  "<int>" → row count (first chunk only).
        // SQL meta directives (`--meta.connectFetchSize: ...`) are mirrored
        // into FuncCall.options under the same bare key, but Datlas may pass
        // them with a `meta.` prefix in some paths, so we check both.
        if let Some(v) = Self::lookup_option(call, "connectFetchSize") {
            if let Some(s) = v.as_str() {
                if let Some(mb) = parse_mb(s) {
                    strat.target_bytes = mb;
                }
                else if let Ok(n) = s.parse::<usize>() {
                    strat.exact_rows = Some(clamp_rows(n));
                }
            }
            else if let Some(n) = v.as_u64() {
                strat.exact_rows = Some(clamp_rows(n as usize));
            }
        }

        if let Some(v) = Self::lookup_option(call, "initConnectFetchSize") {
            if let Some(s) = v.as_str() {
                if let Ok(n) = s.parse::<usize>() {
                    strat.init_exact_rows = Some(clamp_rows(n));
                }
            }
            else if let Some(n) = v.as_u64() {
                strat.init_exact_rows = Some(clamp_rows(n as usize));
            }
        }

        // Per-query override for dictionary encoding of low-cardinality Utf8
        // columns. Values: `off` / `auto` / `<ratio>` (see `DictEncodeConfig::from_str_option`).
        if let Some(v) = Self::lookup_option(call, "dictEncodeStrings") {
            if let Some(s) = v.as_str() {
                strat.dict_encode = DictEncodeConfig::from_str_option(s);
            }
            else if let Some(b) = v.as_bool() {
                strat.dict_encode = if b { Some(DictEncodeConfig::default()) } else { None };
            }
        }

        strat
    }

    /// Look up an option in either FuncCall.options or DataQuery.options, with
    /// or without the `meta.` prefix.
    fn lookup_option<'a>(call: &'a FuncCall, key: &str) -> Option<&'a serde_json::Value> {
        let meta_key = format!("meta.{}", key);
        call.options.get(key)
            .or_else(|| call.options.get(&meta_key))
            .or_else(|| call.func.options.get(key))
            .or_else(|| call.func.options.get(&meta_key))
    }
}

fn parse_mb(s: &str) -> Option<usize> {
    let s = s.trim();
    let n = s.strip_suffix("MB").or_else(|| s.strip_suffix("mb"))?.trim();
    n.parse::<usize>().ok().map(|n| n * 1_000_000)
}

fn clamp_rows(n: usize) -> usize {
    n.clamp(MIN_FETCH_ROWS, MAX_FETCH_ROWS)
}

/// WS /query_socket — streaming query execution with Arrow IPC output.
pub async fn query_socket_upgrade(
    ws: WebSocketUpgrade,
    State(state): State<AppState>,
) -> impl IntoResponse {
    ws.on_upgrade(move |socket| handle_query_socket(socket, state))
}

/// WebSocket protocol state machine (mirrors Java SessionHandler.java):
///
/// 1. Send "CONNECTED"
/// 2. Receive "QUERY <json>" → parse FuncCall, start query execution
/// 3. For each chunk (accumulated RecordBatches → Arrow IPC):
///    a. Send "DATAFRAME PART SIZE: <len> format=arrow-ipc"
///    b. Receive "DATAFRAME PART SIZE RECEIVED" → send binary chunk
///    c. Receive "PART OK"
/// 4. Send "COMPLETED <chunkCount>"
/// 5. Receive "COMPLETED_OK" → send "EOF", close
async fn handle_query_socket(mut socket: WebSocket, state: AppState) {
    if send_text(&mut socket, "CONNECTED").await.is_err() {
        return;
    }

    let func_call = match receive_query(&mut socket).await {
        Ok(call) => call,
        Err(e) => {
            let _ = send_error(&mut socket, &format!("Failed to parse query: {}", e)).await;
            return;
        }
    };

    let data_source = func_call.connection().data_source.clone();
    info!(provider = %data_source, query_id = %func_call.id, "WebSocket query started");

    let provider = match state.registry.get(&data_source) {
        Some(p) => p.clone(),
        None => {
            let _ = send_error(&mut socket, &format!("Unknown data source: {}", data_source)).await;
            return;
        }
    };

    let pool_mgr = state.pool_manager.clone();

    let strategy = ChunkStrategy::from_env_and_call(&func_call);
    info!(
        init_target_bytes = strategy.init_target_bytes,
        target_bytes = strategy.target_bytes,
        max_rows = strategy.max_rows,
        exact_rows = ?strategy.exact_rows,
        init_exact_rows = ?strategy.init_exact_rows,
        dict_encode = ?strategy.dict_encode,
        "Chunk strategy resolved"
    );

    let (tx, rx) = tokio::sync::mpsc::channel(BATCH_CHANNEL_SIZE);
    let query_handle = tokio::task::spawn_blocking(move || {
        provider.execute_query_streaming(&func_call, &pool_mgr, tx)
    });

    let mut chunk_count = 0i32;
    let mut stream_error: Option<String> = None;
    let mut total_bytes: usize = 0;
    let mut t_ws_total: u128 = 0;
    let stream_start = std::time::Instant::now();

    // Pipelined streaming, mirrors Java SessionHandler's CompletableFuture pattern:
    // a separate task accumulates RecordBatches into ready-to-send IPC chunks while
    // the main loop sends the previous chunk over WebSocket. The two run concurrently
    // — ADBC keeps downloading and serializing while WS is busy with the prior chunk.
    let (chunk_tx, mut chunk_rx) = tokio::sync::mpsc::channel::<Vec<u8>>(CHUNK_CHANNEL_SIZE);
    let drain_handle = tokio::spawn(drain_batches_to_chunks(rx, chunk_tx, strategy));

    while let Some(ipc_bytes) = chunk_rx.recv().await {
        total_bytes += ipc_bytes.len();
        let t_ws = std::time::Instant::now();
        if let Err(e) = send_chunk_bytes(&mut socket, ipc_bytes, &mut chunk_count).await {
            error!(error = %e, "Failed to send chunk");
            stream_error = Some(format!("{}", e));
            break;
        }
        t_ws_total += t_ws.elapsed().as_micros();
    }

    let (t_recv_total, t_serialize_total) = drain_handle.await.unwrap_or((0, 0));

    info!(
        total_ms = stream_start.elapsed().as_millis() as u64,
        recv_ms = (t_recv_total / 1000) as u64,
        serialize_ms = (t_serialize_total / 1000) as u64,
        websocket_ms = (t_ws_total / 1000) as u64,
        chunks = chunk_count,
        bytes = total_bytes,
        "Streaming breakdown"
    );

    match query_handle.await {
        Ok(Err(e)) => {
            if stream_error.is_none() {
                let _ = send_error(&mut socket, &format!("{:#}", e)).await;
                return;
            }
        }
        Err(e) => {
            if stream_error.is_none() {
                let _ = send_error(&mut socket, &format!("Internal error: {}", e)).await;
                return;
            }
        }
        Ok(Ok(())) => {}
    }

    if let Some(err) = stream_error {
        let _ = send_error(&mut socket, &err).await;
        return;
    }

    if send_text(&mut socket, &format!("COMPLETED {}", chunk_count)).await.is_err() {
        return;
    }

    match receive_text(&mut socket).await {
        Ok(msg) if msg.starts_with("COMPLETED_OK") => {
            let _ = send_text(&mut socket, "EOF").await;
        }
        Ok(msg) => {
            warn!(msg = %msg, "Expected COMPLETED_OK, got unexpected message");
            let _ = send_text(&mut socket, "EOF").await;
        }
        Err(_) => {}
    }

    info!(chunks = chunk_count, "WebSocket query completed (streaming)");
}

/// Drain RecordBatches from the producer and forward Arrow-IPC-encoded chunks
/// to the WebSocket sender. Flushing mirrors Java `QueryManager`:
///
/// - The first chunk uses the smaller `init_target_bytes` (or
///   `init_exact_rows` if supplied) so the client gets something back fast
///   — this is the TTFR probe.
/// - Subsequent chunks use `target_bytes` (or `exact_rows`), capped by
///   `max_rows` regardless of byte size so tables with tiny rows don't
///   accumulate millions of rows in a single buffer.
///
/// Returns (total recv micros, total serialize micros) for diagnostics.
async fn drain_batches_to_chunks(
    mut rx: tokio::sync::mpsc::Receiver<arrow::record_batch::RecordBatch>,
    chunk_tx: tokio::sync::mpsc::Sender<Vec<u8>>,
    strategy: ChunkStrategy,
) -> (u128, u128) {
    let mut pending: Vec<arrow::record_batch::RecordBatch> = Vec::new();
    let mut pending_rows: usize = 0;
    let mut pending_bytes: usize = 0;
    let mut serialize_total: u128 = 0;
    let mut recv_total: u128 = 0;
    let mut sent_first_chunk = false;

    loop {
        let t = std::time::Instant::now();
        let batch_opt = rx.recv().await;
        recv_total += t.elapsed().as_micros();

        let Some(batch) = batch_opt else { break };
        pending_rows += batch.num_rows();
        pending_bytes += batch.get_array_memory_size();
        pending.push(batch);

        // If the caller set an exact row count (via `connectFetchSize` /
        // `initConnectFetchSize` as an integer), honour it strictly and
        // ignore the byte target — this matches Java `QueryManager`'s
        // `setFetchSize` row-form path and is how we get apples-to-apples
        // benchmark comparisons. Otherwise fall back to the byte target
        // capped by `max_rows`.
        let should_flush = if !sent_first_chunk {
            if let Some(rows) = strategy.init_exact_rows {
                pending_rows >= rows
            }
            else {
                pending_bytes >= strategy.init_target_bytes
                    || pending_rows >= strategy.max_rows
            }
        }
        else {
            if let Some(rows) = strategy.exact_rows {
                pending_rows >= rows
            }
            else {
                pending_bytes >= strategy.target_bytes
                    || pending_rows >= strategy.max_rows
            }
        };

        if should_flush {
            let t = std::time::Instant::now();
            let ipc_bytes = match serialize_pending(&pending, strategy.dict_encode) {
                Ok(b) => b,
                Err(_) => break,
            };
            serialize_total += t.elapsed().as_micros();
            pending.clear();
            pending_rows = 0;
            pending_bytes = 0;
            if chunk_tx.send(ipc_bytes).await.is_err() {
                return (recv_total, serialize_total);
            }
            sent_first_chunk = true;
        }
    }

    if !pending.is_empty() {
        let t = std::time::Instant::now();
        if let Ok(ipc_bytes) = serialize_pending(&pending, strategy.dict_encode) {
            serialize_total += t.elapsed().as_micros();
            let _ = chunk_tx.send(ipc_bytes).await;
        }
    }
    (recv_total, serialize_total)
}

/// Dict-encode and serialise a group of pending batches to an Arrow IPC stream.
///
/// When dict encoding is enabled, the group is first concatenated into a single
/// batch and then encoded once. This is necessary for correctness (each batch in
/// the group must end up with the same schema, otherwise `StreamWriter::write`
/// rejects the second batch as schema-mismatched) and it also gives a single
/// dictionary per chunk instead of one-per-driver-batch, which maximises the
/// compression win. The concat is a single Arrow `concat_batches` call — it
/// memcpys the underlying buffers but avoids per-row work.
fn serialize_pending(
    pending: &[arrow::record_batch::RecordBatch],
    dict_encode: Option<DictEncodeConfig>,
) -> anyhow::Result<Vec<u8>> {
    if pending.is_empty() {
        return Ok(Vec::new());
    }
    if let Some(cfg) = dict_encode {
        let schema = pending[0].schema();
        let merged = arrow::compute::concat_batches(&schema, pending)?;
        let encoded = maybe_encode_batch(&merged, cfg)?;
        crate::arrow_ipc::record_batches_to_ipc(std::slice::from_ref(&encoded))
    }
    else {
        crate::arrow_ipc::record_batches_to_ipc(pending)
    }
}

/// Send pre-serialized IPC bytes through the WebSocket protocol.
async fn send_chunk_bytes(
    socket: &mut WebSocket,
    ipc_bytes: Vec<u8>,
    chunk_count: &mut i32,
) -> anyhow::Result<()> {
    let len = ipc_bytes.len();
    send_text(socket, &format!("DATAFRAME PART SIZE: {} format=arrow-ipc", len)).await?;

    let ack = receive_text(socket).await?;
    if !ack.contains("DATAFRAME PART SIZE RECEIVED") {
        anyhow::bail!("Expected DATAFRAME PART SIZE RECEIVED, got: {}", ack);
    }

    socket.send(Message::Binary(ipc_bytes.into())).await?;

    let ok = receive_text(socket).await?;
    if !ok.contains("PART OK") {
        anyhow::bail!("Expected PART OK, got: {}", ok);
    }

    *chunk_count += 1;
    Ok(())
}

async fn send_text(socket: &mut WebSocket, text: &str) -> anyhow::Result<()> {
    socket.send(Message::Text(text.to_string().into())).await?;
    Ok(())
}

async fn send_error(socket: &mut WebSocket, message: &str) -> anyhow::Result<()> {
    let result = DataQueryRunResult::error(message);
    let json = serde_json::to_string(&result)?;
    send_text(socket, &format!("ERROR: {}", json)).await
}

async fn receive_text(socket: &mut WebSocket) -> anyhow::Result<String> {
    loop {
        match socket.recv().await {
            Some(Ok(Message::Text(text))) => return Ok(text.to_string()),
            Some(Ok(Message::Ping(_))) => continue,
            Some(Ok(Message::Pong(_))) => continue,
            Some(Ok(Message::Close(_))) => anyhow::bail!("WebSocket closed by client"),
            Some(Ok(_)) => continue,
            Some(Err(e)) => anyhow::bail!("WebSocket error: {}", e),
            None => anyhow::bail!("WebSocket closed"),
        }
    }
}

async fn receive_query(socket: &mut WebSocket) -> anyhow::Result<FuncCall> {
    let msg = receive_text(socket).await?;
    let json_str = msg.strip_prefix("QUERY ")
        .ok_or_else(|| anyhow::anyhow!("Expected 'QUERY <json>', got: {}", &msg[..msg.len().min(100)]))?;
    let call: FuncCall = serde_json::from_str(json_str)?;
    Ok(call)
}
