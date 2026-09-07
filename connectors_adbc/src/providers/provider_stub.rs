//! A provider backed by nothing, for tests of registry-driven behaviour — routing,
//! advertising, lookup — that are indifferent to any particular database. Reaching for a
//! real provider there would couple provider-agnostic code to one vendor.

use crate::protocol::{default_engine, DataConnection, DataSourceDescriptor, FuncCall};
use crate::providers::{ConnectionParams, DbProvider};
use crate::query::params::SqlDialect;
use crate::state::{ConnectionPoolManager, PoolDatabase};
use arrow::record_batch::RecordBatch;
use anyhow::Result;
use std::collections::HashMap;

pub struct StubProvider {
    name: String,
}

impl StubProvider {
    pub fn new(name: &str) -> Self {
        StubProvider { name: name.to_string() }
    }
}

impl SqlDialect for StubProvider {}

impl DbProvider for StubProvider {
    fn name(&self) -> &str {
        &self.name
    }

    fn descriptor(&self) -> DataSourceDescriptor {
        DataSourceDescriptor {
            type_tag: "DataSource".to_string(),
            source_type: self.name.clone(),
            engine: default_engine(),
            category: "Database".to_string(),
            description: String::new(),
            comment_start: "--".to_string(),
            query_language: "sql".to_string(),
            can_browse_schema: false,
            default_schema: None,
            name_brackets: "\"".to_string(),
            limit_at_end: true,
            connection_template: Vec::new(),
            credentials_template: Vec::new(),
            aggregations: Vec::new(),
            types_map: HashMap::new(),
        }
    }

    /// Backed by no driver, so nothing to load.
    fn probe(&self) -> Result<()> {
        Ok(())
    }

    fn build_pool_database(&self, _params: &ConnectionParams) -> Result<PoolDatabase> {
        unimplemented!("StubProvider never connects")
    }

    fn connection_params(&self, _conn: &DataConnection) -> Result<ConnectionParams> {
        Ok(ConnectionParams { options: HashMap::new() })
    }

    fn test_connection(&self, _conn: &DataConnection, _pool_mgr: &ConnectionPoolManager) -> Result<()> {
        unimplemented!("StubProvider never connects")
    }

    fn execute_query_streaming(
        &self,
        _call: &FuncCall,
        _pool_mgr: &ConnectionPoolManager,
        _tx: tokio::sync::mpsc::Sender<RecordBatch>,
    ) -> Result<()> {
        unimplemented!("StubProvider never queries")
    }
}
