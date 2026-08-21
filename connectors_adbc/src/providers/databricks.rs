use crate::protocol::{default_engine, DataConnection, DataSourceDescriptor, FuncCall, Property};
use crate::providers::{Advertised, ConnectionParams, DbProvider};
use crate::query::params::{SqlDialect, prepare_query};
use crate::state::{ConnectionPoolManager, PoolDatabase};
use adbc_core::options::{OptionDatabase, OptionValue};
use adbc_core::{Driver as DriverTrait, Optionable};
use arrow::record_batch::RecordBatch;
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use tracing::info;

// ADBC option keys (from databricks-adbc driver docs).
const OPT_URI: &str = "uri";
const OPT_HTTP_PATH: &str = "databricks.http_path";
const OPT_AUTH_TYPE: &str = "databricks.auth.type";
const OPT_ACCESS_TOKEN: &str = "databricks.access_token";
const OPT_CATALOG: &str = "databricks.catalog";
const OPT_CF_MAX_CHUNKS: &str = "databricks.cloudfetch.max_chunks_in_memory";
const OPT_CF_PREFETCH: &str = "databricks.cloudfetch.link_prefetch_window";
const AUTH_ACCESS_TOKEN: &str = "access_token";

pub struct DatabricksProvider;

impl DatabricksProvider {
    pub fn new() -> Self {
        DatabricksProvider
    }

    /// Parse JDBC URL: jdbc:databricks://host:port/catalog;...;httpPath=...
    pub fn parse_jdbc_url(url: &str) -> Result<JdbcUrlParts> {
        let url = url.trim();
        let stripped = url.strip_prefix("jdbc:databricks://")
            .context("JDBC URL must start with jdbc:databricks://")?;

        // Split host:port/catalog from parameters
        let (host_part, params_str) = match stripped.find(';') {
            Some(i) => (&stripped[..i], &stripped[i + 1..]),
            None => (stripped, ""),
        };

        // Parse host:port/catalog
        let (host_port, catalog) = match host_part.find('/') {
            Some(i) => (&host_part[..i], &host_part[i + 1..]),
            None => (host_part, "default"),
        };

        let (host, port) = match host_port.rfind(':') {
            Some(i) => {
                let h = &host_port[..i];
                let p = host_port[i + 1..].parse::<u16>().unwrap_or(443);
                (h.to_string(), p)
            }
            None => (host_port.to_string(), 443),
        };

        // Parse semicolon-separated parameters
        let mut http_path = String::new();
        for param in params_str.split(';') {
            let param = param.trim();
            if param.is_empty() {
                continue;
            }
            if let Some(i) = param.find('=') {
                let key = param[..i].trim().to_lowercase();
                let val = param[i + 1..].trim();
                if key == "httppath" {
                    http_path = val.to_string();
                }
            }
        }

        Ok(JdbcUrlParts {
            host,
            port,
            catalog: catalog.to_string(),
            http_path,
        })
    }

    /// Extract connection params from a DataConnection.
    /// Supports both Datagrok platform format (workspaceURL, token) and direct params (server, password).
    fn extract_params(conn: &DataConnection) -> Result<DatabricksConnParams> {
        // Try workspaceURL first (Datagrok platform format), then server/host (direct)
        if let Some(host) = conn.get_str("workspaceURL")
            .or_else(|| conn.get_str("server"))
            .or_else(|| conn.get_str("host"))
        {
            let port = conn.get_i64("port").unwrap_or(443) as u16;
            let http_path = conn.get_str("httpPath").unwrap_or("").to_string();
            let catalog = conn.get_str("db").unwrap_or("").to_string();
            // Token: try credential "token" (platform), then "password", then "#token" (federated)
            let token = conn.get_credential("token")
                .or_else(|| conn.get_credential("password"))
                .or_else(|| conn.get_credential("#token"))
                .unwrap_or("")
                .to_string();

            return Ok(DatabricksConnParams { host: host.to_string(), port, http_path, catalog, token });
        }

        // Fall back to JDBC URL in connectionString or connString parameter
        let jdbc_url = conn.connection_string.as_deref()
            .or_else(|| conn.get_str("connString"))
            .or_else(|| conn.get_str("connectionString"))
            .context("No workspaceURL, server, or JDBC URL provided for Databricks connection")?;

        let parts = Self::parse_jdbc_url(jdbc_url)?;
        let token = conn.get_credential("token")
            .or_else(|| conn.get_credential("password"))
            .unwrap_or("")
            .to_string();

        Ok(DatabricksConnParams {
            host: parts.host,
            port: parts.port,
            http_path: parts.http_path,
            catalog: parts.catalog,
            token,
        })
    }
}

#[derive(Debug)]
pub struct JdbcUrlParts {
    pub host: String,
    pub port: u16,
    pub catalog: String,
    pub http_path: String,
}

#[derive(Debug)]
struct DatabricksConnParams {
    host: String,
    port: u16,
    http_path: String,
    catalog: String,
    token: String,
}

/// Not served: the Java `Databricks` connector handles this type today. Flip to `true` to
/// advertise `DatabricksADBC` alongside it — the two are distinct types, so they coexist.
impl Advertised for DatabricksProvider {
    const ADVERTISED: bool = false;
}

/// Vendor-neutral search-pattern SQL is accepted as-is. Java's `DatabricksProvider` does not
/// override `getRegexQuery` either, so `regex` search patterns are rejected on both.
impl SqlDialect for DatabricksProvider {}

impl DbProvider for DatabricksProvider {
    fn name(&self) -> &str {
        "DatabricksADBC"
    }

    /// The Databricks driver is a Rust crate linked into this binary, not a `dlopen`ed
    /// library, so there is nothing that can fail to load at startup.
    fn probe(&self) -> Result<()> {
        Ok(())
    }

    fn build_pool_database(&self, params: &ConnectionParams) -> Result<PoolDatabase> {
        let mut driver = databricks_adbc::Driver::new();
        let mut db = DriverTrait::new_database(&mut driver)
            .context("Failed to create Databricks ADBC database")?;
        for (key_name, value) in &params.options {
            let opt = if key_name == "uri" {
                OptionDatabase::Uri
            }
            else {
                OptionDatabase::Other(key_name.clone().into())
            };
            db.set_option(opt, OptionValue::String(value.clone()))
                .with_context(|| format!("Failed to set Databricks option: {}", key_name))?;
        }
        Ok(Box::new(db))
    }

    fn descriptor(&self) -> DataSourceDescriptor {
        DataSourceDescriptor {
            type_tag: "DataSource".to_string(),
            source_type: "DatabricksADBC".to_string(),
            engine: default_engine(),
            category: "Database".to_string(),
            description: "Databricks SQL via ADBC (Arrow IPC)".to_string(),
            comment_start: "--".to_string(),
            query_language: "sql".to_string(),
            can_browse_schema: false,
            default_schema: None,
            name_brackets: "`".to_string(),
            limit_at_end: true,
            connection_template: vec![
                Property::string("server").with_description("Databricks workspace hostname"),
                Property::int("port").with_description("Port (default: 443)"),
                Property::string("httpPath").with_description("HTTP path to SQL warehouse"),
                Property::string("db").with_description("Default catalog"),
            ],
            credentials_template: vec![
                Property::string("password").as_password().with_description("Personal Access Token (PAT)"),
            ],
            aggregations: vec![],
            types_map: HashMap::new(),
        }
    }

    fn connection_params(&self, conn: &DataConnection) -> Result<ConnectionParams> {
        let p = Self::extract_params(conn)?;
        let mut options = HashMap::new();
        options.insert(OPT_URI.to_string(), format!("https://{}:{}", p.host, p.port));
        options.insert(OPT_HTTP_PATH.to_string(), p.http_path);
        options.insert(OPT_AUTH_TYPE.to_string(), AUTH_ACCESS_TOKEN.to_string());
        options.insert(OPT_ACCESS_TOKEN.to_string(), p.token);
        if !p.catalog.is_empty() {
            options.insert(OPT_CATALOG.to_string(), p.catalog);
        }
        // CloudFetch tuning: matches Java JDBC's cloudFetchThreadPoolSize=16 for S3 throughput.
        options.insert(OPT_CF_MAX_CHUNKS.to_string(), "16".to_string());
        options.insert(OPT_CF_PREFETCH.to_string(), "128".to_string());
        Ok(ConnectionParams { options })
    }

    fn test_connection(&self, conn: &DataConnection, pool_mgr: &ConnectionPoolManager) -> Result<()> {
        let params = Self::extract_params(conn)?;
        info!(host = %params.host, "Testing Databricks connection");

        let conn_params = self.connection_params(conn)?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, conn, &conn_params)?;

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query("SELECT 1")?;
        let batches = stmt.execute_to_batches()?;
        if batches.is_empty() || batches[0].num_rows() < 1 {
            bail!("SELECT 1 returned no rows");
        }
        info!("Databricks connection test successful");
        Ok(())
    }

    fn execute_query_streaming(
        &self,
        call: &FuncCall,
        pool_mgr: &ConnectionPoolManager,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
    ) -> Result<()> {
        let prepared = prepare_query(call, self)?;
        info!(sql = %prepared.sql, "Executing Databricks query (streaming)");

        let conn_params = self.connection_params(call.connection())?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, call.connection(), &conn_params)?;

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query(&prepared.sql)?;
        stmt.execute_streaming(tx)?;
        info!("Databricks streaming query completed");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_jdbc_url_standard() {
        let url = "jdbc:databricks://dbc-52240d8b-a70a.cloud.databricks.com:443/default;transportMode=http;ssl=1;AuthMech=3;httpPath=/sql/1.0/warehouses/769f690b943edd7f";
        let parts = DatabricksProvider::parse_jdbc_url(url).unwrap();
        assert_eq!(parts.host, "dbc-52240d8b-a70a.cloud.databricks.com");
        assert_eq!(parts.port, 443);
        assert_eq!(parts.catalog, "default");
        assert_eq!(parts.http_path, "/sql/1.0/warehouses/769f690b943edd7f");
    }

    #[test]
    fn test_parse_jdbc_url_no_port() {
        let url = "jdbc:databricks://myhost.databricks.com/mycatalog;httpPath=/path";
        let parts = DatabricksProvider::parse_jdbc_url(url).unwrap();
        assert_eq!(parts.host, "myhost.databricks.com");
        assert_eq!(parts.port, 443);
        assert_eq!(parts.catalog, "mycatalog");
        assert_eq!(parts.http_path, "/path");
    }

    #[test]
    fn test_parse_jdbc_url_no_catalog() {
        let url = "jdbc:databricks://myhost.databricks.com:443;httpPath=/p";
        let parts = DatabricksProvider::parse_jdbc_url(url).unwrap();
        assert_eq!(parts.host, "myhost.databricks.com");
        assert_eq!(parts.port, 443);
        assert_eq!(parts.catalog, "default");
    }

    #[test]
    fn test_parse_jdbc_url_invalid() {
        let result = DatabricksProvider::parse_jdbc_url("jdbc:postgresql://localhost/db");
        assert!(result.is_err());
    }
}
