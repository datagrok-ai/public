use crate::protocol::{DataConnection, DataSourceDescriptor, FuncCall};
use crate::query::params::SqlDialect;
use crate::state::{ConnectionPoolManager, PoolDatabase};
use anyhow::Result;
use arrow::record_batch::RecordBatch;
use std::collections::HashMap;

/// Database provider trait — intentionally sync.
/// All methods block. Callers wrap in tokio::task::spawn_blocking().
///
/// A provider is its own [`SqlDialect`]: search-pattern rendering is expressed by
/// overriding trait methods in the provider's own file, so nothing outside
/// `providers/` learns about a specific database.
pub trait DbProvider: SqlDialect + Send + Sync {
    fn name(&self) -> &str;
    fn descriptor(&self) -> DataSourceDescriptor;

    /// Build the driver-specific ADBC database for a new connection pool.
    /// Each provider owns its driver loading and connection-option mapping.
    fn build_pool_database(&self, params: &ConnectionParams) -> Result<PoolDatabase>;

    /// Verify this provider can actually execute queries: load its ADBC driver and drop it.
    /// Run once at startup, before registration — see [`build_registry`](super::build_registry).
    /// Only advertised providers are probed.
    ///
    /// Datlas routes every query for an advertised type to whichever endpoint claims it, so a
    /// provider that is advertised but whose driver won't load silently swallows that type
    /// rather than leaving it to the Java default. Driver loading otherwise happens lazily in
    /// [`DbProvider::build_pool_database`], which would surface the failure only at the first
    /// query — long after `/conn` promised we could serve it.
    ///
    /// Providers whose driver is linked in at compile time are always available and return
    /// `Ok(())`; only the dynamically-loaded ones can fail here.
    fn probe(&self) -> Result<()>;

    /// Build ADBC connection parameters for creating a new pool.
    /// Returns (driver_init_func, connection_options).
    fn connection_params(&self, conn: &DataConnection) -> Result<ConnectionParams>;

    /// What partitions this provider's connection pool, within the provider's own namespace.
    ///
    /// Connections sharing a key share ADBC connections, and so the credentials behind them —
    /// a provider whose connections the default identity cannot tell apart must override this.
    fn pool_identity(&self, conn: &DataConnection) -> String {
        conn.connection_identity()
    }

    fn test_connection(&self, conn: &DataConnection, pool_mgr: &ConnectionPoolManager) -> Result<()>;

    fn schemas_sql(&self) -> Option<String> {
        None
    }

    fn schema_sql(&self, _schema: Option<&str>, _table: Option<&str>, _include_key_info: bool) -> Option<String> {
        None
    }

    /// Stream query results through a channel. Each RecordBatch is sent as it arrives.
    /// This is the only method a provider implements — avoids collecting everything in memory.
    fn execute_query_streaming(
        &self,
        call: &FuncCall,
        pool_mgr: &ConnectionPoolManager,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
    ) -> Result<()>;

    /// Collect the whole result set by draining `execute_query_streaming`.
    /// Blocking, and must not be called from an async context (the streaming side
    /// uses `blocking_send`); the server always streams instead.
    fn execute_query(&self, call: &FuncCall, pool_mgr: &ConnectionPoolManager) -> Result<Vec<RecordBatch>> {
        let (tx, mut rx) = tokio::sync::mpsc::channel(COLLECT_CHANNEL_CAPACITY);
        std::thread::scope(|scope| {
            let producer = scope.spawn(|| self.execute_query_streaming(call, pool_mgr, tx));
            let mut batches = Vec::new();
            while let Some(batch) = rx.blocking_recv() {
                batches.push(batch);
            }
            producer.join().map_err(|_| anyhow::anyhow!("query thread panicked"))??;
            Ok(batches)
        })
    }
}

/// Backpressure window while collecting a streamed result set in memory.
const COLLECT_CHANNEL_CAPACITY: usize = 32;

/// Whether [`build_registry`](super::build_registry) serves a provider.
///
/// Deliberately not a `DbProvider` method: this is a fact about the provider type, not
/// something a live instance answers, and `DbProvider` is used as `dyn DbProvider`
/// throughout — a trait carrying associated consts is not object-safe, so it cannot live
/// there. `register_probed` takes it as a bound instead, where the concrete type is known.
pub trait Advertised {
    /// `false` keeps the provider compiled, tested and ready to re-enable, but out of the
    /// registry — so `GET /conn` never claims the type and datlas leaves it to the Java
    /// default. Its driver is not probed either: an unserved provider must not be able to
    /// stop the service from starting.
    const ADVERTISED: bool;
}

/// Parameters needed to establish an ADBC connection.
/// Used by the pool manager to create new connections.
pub struct ConnectionParams {
    /// Key-value pairs for ADBC database/connection options.
    pub options: HashMap<String, String>,
}
