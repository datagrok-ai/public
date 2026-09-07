mod arrow_ipc;
mod config;
mod dict_encode;
mod protocol;
mod providers;
mod query;
mod routes;
mod state;

use crate::config::Config;
use crate::providers::build_registry;
use crate::state::{AppState, ConnectionPoolManager};
use axum::routing::{get, post};
use axum::Router;
use std::io::{Read, Write};
use std::sync::Arc;
use std::time::Duration;
use tracing::info;

/// Prepend `<crate-root>/lib` to the OS dynamic-library search path so that
/// `ManagedDriver::load_dynamic_from_name("adbc_driver_snowflake", ...)` finds
/// the fetched `adbc_driver_snowflake.{dll,so,dylib}` without requiring a
/// system-wide install.
///
/// Safety: this must run before any tokio worker or adbc driver-load happens.
/// We invoke it at the very top of `main()`, before `tracing_subscriber::fmt`
/// and before the tokio runtime starts spawning threads.
fn prepend_driver_search_path() {
    let lib = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("lib");
    if !lib.exists() {
        return;
    }
    #[cfg(target_os = "windows")]
    let var = "PATH";
    #[cfg(target_os = "linux")]
    let var = "LD_LIBRARY_PATH";
    #[cfg(target_os = "macos")]
    let var = "DYLD_LIBRARY_PATH";
    let old = std::env::var(var).unwrap_or_default();
    let sep = if cfg!(windows) { ';' } else { ':' };
    let new_value = if old.is_empty() {
        lib.display().to_string()
    }
    else {
        format!("{}{}{}", lib.display(), sep, old)
    };
    std::env::set_var(var, new_value);
}

/// Probes the already-running server's `/health` and exits with its verdict. The Docker
/// HEALTHCHECK calls this instead of curl, which the distroless runtime image does not carry.
/// Hand-rolled over a raw socket so the server binary gains no HTTP-client dependency.
fn health_probe(port: u16) -> ! {
    let probe = || -> std::io::Result<bool> {
        let mut socket = std::net::TcpStream::connect(("127.0.0.1", port))?;
        socket.set_read_timeout(Some(Duration::from_secs(3)))?;
        socket.write_all(b"GET /health HTTP/1.1\r\nHost: localhost\r\nConnection: close\r\n\r\n")?;
        let mut status = [0u8; 12];
        socket.read_exact(&mut status)?;
        Ok(&status == b"HTTP/1.1 200")
    };
    std::process::exit(if probe().unwrap_or(false) { 0 } else { 1 });
}

#[tokio::main]
async fn main() {
    if std::env::args().any(|a| a == "--health") {
        health_probe(Config::from_env().port);
    }

    prepend_driver_search_path();

    // Initialize tracing
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "grok_connect_adbc=info,tower_http=info".into()),
        )
        .init();

    let config = Config::from_env();
    config.log_config();

    // Probes each provider's driver; exits non-zero if one won't load.
    let registry = build_registry();

    // Build connection pool manager
    let pool_manager = Arc::new(ConnectionPoolManager::new(
        config.pool_max_size,
        config.pool_idle_timeout_secs,
    ));

    // Spawn background pool eviction task
    let eviction_pool_mgr = pool_manager.clone();
    let eviction_interval = config.pool_eviction_interval_secs;
    tokio::spawn(async move {
        let mut interval = tokio::time::interval(Duration::from_secs(eviction_interval));
        loop {
            interval.tick().await;
            eviction_pool_mgr.evict_idle_pools();
        }
    });

    let state = AppState {
        registry: Arc::new(registry),
        pool_manager,
        config: Arc::new(config.clone()),
    };

    // Build router
    let app = Router::new()
        .route("/health", get(routes::health::health))
        .route("/info", get(routes::health::info))
        .route("/conn", get(routes::conn::list_providers))
        .route("/test", post(routes::test_connection::test_connection))
        .route("/query_socket", get(routes::query_socket::query_socket_upgrade))
        .route("/query", post(routes::query::query_stub))
        .route("/schema", post(routes::schema::schema))
        .route("/schemas", post(routes::schema::schemas))
        .route("/cancel", post(routes::cancel::cancel_stub))
        .with_state(state);

    let addr = format!("0.0.0.0:{}", config.port);
    info!(addr = %addr, "Starting grok_connect_adbc server");

    let listener = tokio::net::TcpListener::bind(&addr).await.unwrap();
    axum::serve(listener, app).await.unwrap();
}
