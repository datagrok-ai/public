use crate::config::Config;
use crate::providers::{ConnectionParams, DbProvider, ProviderRegistry};
use crate::protocol::DataConnection;
use arrow::record_batch::RecordBatch;
use anyhow::{Result, Context};
use std::collections::HashMap;
use std::sync::{Arc, RwLock, Mutex};
use std::time::Duration;
use tracing::info;

use adbc_core::{Connection as ConnTrait, Database as DbTrait, Statement as StmtTrait};

/// Shared application state, passed to all route handlers.
#[derive(Clone)]
pub struct AppState {
    pub registry: Arc<ProviderRegistry>,
    pub pool_manager: Arc<ConnectionPoolManager>,
    pub config: Arc<Config>,
}

// Object-safe wrappers that erase each ADBC driver's concrete
// `Database`/`Connection`/`Statement` types — and their associated types, which
// make `adbc_core`'s own traits non-`dyn`-able — behind trait objects. The pool
// holds `Box<dyn Pooled*>` without ever naming a driver; the blanket impls below
// cover every `adbc_core` driver, so a provider just `Box::new(db)`s its database.

pub trait PooledDatabase {
    fn new_connection(&self) -> Result<Box<dyn PooledConnection>>;
}

pub trait PooledConnection {
    fn new_statement(&mut self) -> Result<Box<dyn PooledStatement>>;
}

pub trait PooledStatement {
    fn set_sql_query(&mut self, query: &str) -> Result<()>;
    /// Execute a statement that returns no result set (e.g. `SET ...`).
    fn execute_update(&mut self) -> Result<()>;
    fn execute_to_batches(&mut self) -> Result<Vec<RecordBatch>>;
    /// Stream results through `tx`, rewriting each batch with `transform`.
    fn execute_streaming(
        &mut self,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
        transform: &dyn Fn(RecordBatch) -> Result<RecordBatch>,
    ) -> Result<()>;
}

impl<D> PooledDatabase for D
where
    D: DbTrait,
    D::ConnectionType: 'static,
    <D::ConnectionType as ConnTrait>::StatementType: 'static,
{
    fn new_connection(&self) -> Result<Box<dyn PooledConnection>> {
        Ok(Box::new(DbTrait::new_connection(self).context("Failed to create ADBC connection")?))
    }
}

impl<C> PooledConnection for C
where
    C: ConnTrait,
    C::StatementType: 'static,
{
    fn new_statement(&mut self) -> Result<Box<dyn PooledStatement>> {
        Ok(Box::new(ConnTrait::new_statement(self).context("Failed to create ADBC statement")?))
    }
}

impl<S> PooledStatement for S
where S: StmtTrait
{
    fn set_sql_query(&mut self, query: &str) -> Result<()> {
        StmtTrait::set_sql_query(self, query)?;
        Ok(())
    }

    fn execute_update(&mut self) -> Result<()> {
        StmtTrait::execute_update(self)?;
        Ok(())
    }

    fn execute_to_batches(&mut self) -> Result<Vec<RecordBatch>> {
        let mut batches = Vec::new();
        for batch_result in StmtTrait::execute(self)? {
            batches.push(batch_result.map_err(|e| anyhow::anyhow!("Arrow error: {}", e))?);
        }
        Ok(batches)
    }

    fn execute_streaming(
        &mut self,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
        transform: &dyn Fn(RecordBatch) -> Result<RecordBatch>,
    ) -> Result<()> {
        pump(StmtTrait::execute(self)?, tx, transform)
    }
}

/// Driver-erased ADBC handles held by the pool (see the trait docs above).
pub type PoolDatabase = Box<dyn PooledDatabase>;
type PoolConnection = Box<dyn PooledConnection>;
type PoolStatement = Box<dyn PooledStatement>;

/// A live ADBC connection handed out from the pool.
/// On drop, the underlying connection is returned to the pool for reuse.
pub struct AdbcConnection {
    inner: Option<PoolConnection>,
    pool_entry: Option<Arc<PoolEntry>>,
}

impl AdbcConnection {
    pub fn new_statement(&mut self) -> Result<AdbcStatementWrapper> {
        let conn = self.inner.as_mut().expect("connection already returned to pool");
        Ok(AdbcStatementWrapper { inner: conn.new_statement()? })
    }
}

impl Drop for AdbcConnection {
    fn drop(&mut self) {
        if let (Some(conn), Some(entry)) = (self.inner.take(), self.pool_entry.take()) {
            let mut idle = entry.idle_connections.lock().unwrap();
            if idle.len() < entry.max_size as usize {
                idle.push(conn);
            }
            entry.active_count.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
        }
    }
}

pub struct AdbcStatementWrapper {
    inner: PoolStatement,
}

impl AdbcStatementWrapper {
    pub fn set_sql_query(&mut self, query: &str) -> Result<()> {
        self.inner.set_sql_query(query)
    }

    /// Execute a statement that returns no result set (e.g. `SET ...`).
    pub fn execute_update(&mut self) -> Result<()> {
        self.inner.execute_update()
    }

    pub fn execute_to_batches(&mut self) -> Result<Vec<RecordBatch>> {
        self.inner.execute_to_batches()
    }

    /// Stream results through a channel. Each RecordBatch is sent as it arrives.
    pub fn execute_streaming(
        &mut self,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
    ) -> Result<()> {
        self.execute_streaming_with(tx, Ok)
    }

    /// Stream results through a channel, rewriting each RecordBatch on the way out.
    /// Providers use this to repair types their driver reports imprecisely.
    pub fn execute_streaming_with(
        &mut self,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
        transform: impl Fn(RecordBatch) -> Result<RecordBatch>,
    ) -> Result<()> {
        self.inner.execute_streaming(tx, &transform)
    }
}

/// Drain an ADBC reader into `tx`, applying `transform` and skipping empty batches.
fn pump<E: std::fmt::Display>(
    reader: impl Iterator<Item = std::result::Result<RecordBatch, E>>,
    tx: tokio::sync::mpsc::Sender<RecordBatch>,
    transform: impl Fn(RecordBatch) -> Result<RecordBatch>,
) -> Result<()> {
    for batch_result in reader {
        let batch = batch_result.map_err(|e| anyhow::anyhow!("Arrow error: {}", e))?;
        if batch.num_rows() == 0 {
            continue;
        }
        if tx.blocking_send(transform(batch)?).is_err() {
            break;
        }
    }
    Ok(())
}

// AdbcConnection is Send — the drivers handle thread-safety internally and we
// always take &mut self on the wrappers.
unsafe impl Send for AdbcConnection {}
unsafe impl Sync for AdbcConnection {}

type PoolKey = String;

/// Connection pool manager with real connection reuse.
/// Idle connections are returned to the pool on drop and reused by subsequent requests.
pub struct ConnectionPoolManager {
    pools: RwLock<HashMap<PoolKey, Arc<PoolEntry>>>,
    max_size_per_pool: u32,
    idle_timeout: Duration,
}

struct PoolEntry {
    database: RwLock<PoolDatabase>,
    idle_connections: Mutex<Vec<PoolConnection>>,
    active_count: std::sync::atomic::AtomicUsize,
    max_size: u32,
}

unsafe impl Send for PoolEntry {}
unsafe impl Sync for PoolEntry {}

impl ConnectionPoolManager {
    pub fn new(max_size: u32, idle_timeout_secs: u64) -> Self {
        ConnectionPoolManager {
            pools: RwLock::new(HashMap::new()),
            max_size_per_pool: max_size,
            idle_timeout: Duration::from_secs(idle_timeout_secs),
        }
    }

    pub fn get_or_create_connection(
        &self,
        provider: &dyn DbProvider,
        conn: &DataConnection,
        params: &ConnectionParams,
    ) -> Result<AdbcConnection> {
        let key = self.make_key(provider, conn);

        // Fast path: pool already exists
        {
            let pools = self.pools.read().unwrap();
            if let Some(entry) = pools.get(&key) {
                return Self::get_from_entry(entry.clone());
            }
        }

        // Slow path: create new pool
        let mut pools = self.pools.write().unwrap();
        if let Some(entry) = pools.get(&key) {
            return Self::get_from_entry(entry.clone());
        }

        let database = provider.build_pool_database(params)?;

        let entry = Arc::new(PoolEntry {
            database: RwLock::new(database),
            idle_connections: Mutex::new(Vec::new()),
            active_count: std::sync::atomic::AtomicUsize::new(0),
            max_size: self.max_size_per_pool,
        });
        pools.insert(key.clone(), entry.clone());

        info!(key = %key, "Created new connection pool");
        Self::get_from_entry(entry)
    }

    fn get_from_entry(entry: Arc<PoolEntry>) -> Result<AdbcConnection> {
        // Reuse an idle connection if available
        let reused = {
            let mut idle = entry.idle_connections.lock().unwrap();
            idle.pop()
        };

        if let Some(conn) = reused {
            entry.active_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return Ok(AdbcConnection {
                inner: Some(conn),
                pool_entry: Some(entry),
            });
        }

        // No idle connection — create a new one from the pool's database.
        let connection = {
            let db = entry.database.read().unwrap();
            db.new_connection()?
        };
        entry.active_count.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        Ok(AdbcConnection {
            inner: Some(connection),
            pool_entry: Some(entry),
        })
    }

    fn make_key(&self, provider: &dyn DbProvider, conn: &DataConnection) -> String {
        format!("{}::{}", provider.name(), provider.pool_identity(conn))
    }

    pub fn evict_idle_pools(&self) {
        let mut pools = self.pools.write().unwrap();
        pools.retain(|key, entry| {
            let idle = entry.idle_connections.lock().unwrap();
            let active = entry.active_count.load(std::sync::atomic::Ordering::Relaxed);
            if idle.is_empty() && active == 0 {
                info!(key = %key, "Evicting idle pool");
                false
            }
            else {
                true
            }
        });
    }
}
