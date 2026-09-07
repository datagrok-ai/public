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

// ADBC option keys for the Snowflake driver.
// See: apache/arrow-adbc go/adbc/driver/snowflake docs.
const OPT_ACCOUNT: &str = "adbc.snowflake.sql.account";
const OPT_DB: &str = "adbc.snowflake.sql.db";
const OPT_WAREHOUSE: &str = "adbc.snowflake.sql.warehouse";
const OPT_SCHEMA: &str = "adbc.snowflake.sql.schema";
const OPT_ROLE: &str = "adbc.snowflake.sql.role";
// When false, the Snowflake ADBC driver returns NUMBER(_,0) as Int64 and
// NUMBER(_,>0) as Float64 instead of Decimal128. The Datagrok Arrow JS decoder
// (`public/packages/Arrow/src/api/api.ts::fromFeather`) doesn't handle Decimal
// — `vector.toArray()` returns a typed array with length != rowCount, which
// downstream trips DataFrame.fromColumns with "row counts differ". Keeping
// high_precision off is also a better match for Datagrok's column types
// (Float64 / Int64) and avoids loss-of-precision surprises down the line.
const OPT_USE_HIGH_PRECISION: &str = "adbc.snowflake.sql.client_option.use_high_precision";
// adbc_core exposes USERNAME/PASSWORD as typed enum variants; using the bare
// string keys lets us round-trip through the `ConnectionParams` HashMap and rely
// on state.rs to translate them back into `OptionDatabase::Username`/`Password`.
const OPT_USERNAME: &str = "username";
const OPT_PASSWORD: &str = "password";

pub struct SnowflakeProvider;

impl SnowflakeProvider {
    pub fn new() -> Self {
        SnowflakeProvider
    }

    /// The driver is not compiled in — `try_load` dlopens `adbc_driver_snowflake` from the
    /// OS library search path (see `prepend_driver_search_path` in `main.rs`).
    fn load_driver() -> Result<adbc_snowflake::Driver> {
        adbc_snowflake::Driver::try_load().context(
            "Failed to load adbc_driver_snowflake. Run \
             'scripts/fetch_adbc_driver.{sh,ps1} snowflake' from the grok_connect_adbc crate \
             root, or install libadbc-driver-snowflake system-wide."
        )
    }

    /// Build the Snowflake account string from a DataConnection.
    /// Prefers the platform `accountLocator`/`region`/`cloud` triple (mirrors the
    /// Java `SnowflakeDataProvider.buildAccount` logic), falling back to a direct
    /// `server` field for the Rust integration tests.
    fn extract_params(conn: &DataConnection) -> Result<SnowflakeConnParams> {
        let account = if let Some(locator) = conn.get_str("accountLocator") {
            let mut s = locator.to_string();
            if let Some(r) = conn.get_str("region") {
                if !r.is_empty() {
                    s.push('.');
                    s.push_str(r);
                }
            }
            if let Some(c) = conn.get_str("cloud") {
                if !c.is_empty() {
                    s.push('.');
                    s.push_str(c);
                }
            }
            s
        }
        else {
            conn.get_str("server")
                .context("Snowflake: accountLocator or server is required")?
                .to_string()
        };

        let db = conn.get_str("db").unwrap_or("").to_string();
        let warehouse = conn.get_str("warehouse").unwrap_or("").to_string();
        let schema = conn.get_str("schema").unwrap_or("").to_string();
        let role = conn.get_str("role").unwrap_or("").to_string();

        let login = conn.get_credential("login")
            .context("Snowflake: login is required")?
            .to_string();
        let password = conn.get_credential("password")
            .context("Snowflake: password is required")?
            .to_string();

        Ok(SnowflakeConnParams { account, db, warehouse, schema, role, login, password })
    }
}

#[derive(Debug)]
struct SnowflakeConnParams {
    account: String,
    db: String,
    warehouse: String,
    schema: String,
    role: String,
    login: String,
    password: String,
}

/// Not served: the Java `Snowflake` connector handles this type today. Flip to `true` to
/// advertise `SnowflakeADBC` alongside it — the two are distinct types, so they coexist.
/// Requires `adbc_driver_snowflake` to be present, or the startup probe fails.
impl Advertised for SnowflakeProvider {
    const ADVERTISED: bool = false;
}

/// Snowflake spells regex matching with the `REGEXP` infix operator (mirrors Java's
/// `SnowflakeDataProvider.getRegexQuery`); the rest is the vendor-neutral default.
impl SqlDialect for SnowflakeProvider {
    fn regex_match(&self, column: &str, expression: &str) -> Result<String> {
        Ok(format!("{} REGEXP '{}'", column, expression.replace('\'', "''")))
    }
}

impl DbProvider for SnowflakeProvider {
    fn name(&self) -> &str {
        "SnowflakeADBC"
    }

    fn probe(&self) -> Result<()> {
        Self::load_driver().map(|_| ())
    }

    fn build_pool_database(&self, params: &ConnectionParams) -> Result<PoolDatabase> {
        let mut driver = Self::load_driver()?;
        let mut db = DriverTrait::new_database(&mut driver)
            .context("Failed to create Snowflake ADBC database")?;
        for (key_name, value) in &params.options {
            let opt = if key_name == "uri" {
                OptionDatabase::Uri
            }
            else if key_name == "username" {
                OptionDatabase::Username
            }
            else if key_name == "password" {
                OptionDatabase::Password
            }
            else {
                OptionDatabase::Other(key_name.clone().into())
            };
            db.set_option(opt, OptionValue::String(value.clone()))
                .with_context(|| format!("Failed to set Snowflake option: {}", key_name))?;
        }
        Ok(Box::new(db))
    }

    fn descriptor(&self) -> DataSourceDescriptor {
        let mut types_map = HashMap::new();
        // Mirrors Java SnowflakeDataProvider.typesMap (line 178-191).
        types_map.insert("number".to_string(), "double".to_string());
        types_map.insert("float".to_string(), "double".to_string());
        types_map.insert("text".to_string(), "string".to_string());
        types_map.insert("boolean".to_string(), "bool".to_string());
        types_map.insert("date".to_string(), "datetime".to_string());
        types_map.insert("time".to_string(), "datetime".to_string());
        types_map.insert("#timestamp.*".to_string(), "datetime".to_string());
        types_map.insert("array".to_string(), "list".to_string());
        types_map.insert("object".to_string(), "object".to_string());
        types_map.insert("variant".to_string(), "object".to_string());
        types_map.insert("geography".to_string(), "object".to_string());
        types_map.insert("binary".to_string(), "blob".to_string());

        DataSourceDescriptor {
            type_tag: "DataSource".to_string(),
            source_type: "SnowflakeADBC".to_string(),
            engine: default_engine(),
            category: "Database".to_string(),
            description: "Snowflake via ADBC (Arrow IPC)".to_string(),
            comment_start: "--".to_string(),
            query_language: "sql".to_string(),
            can_browse_schema: false,
            default_schema: Some("PUBLIC".to_string()),
            name_brackets: "\"".to_string(),
            limit_at_end: true,
            connection_template: vec![
                Property::string("accountLocator").with_description("Account locator (e.g., WM00888)"),
                Property::string("region").with_description("Region (e.g., us-east-2)"),
                Property::string("cloud").with_description("Cloud provider (aws / azure / gcp)"),
                Property::string("db").with_description("Database name"),
                Property::string("warehouse").with_description("Warehouse name"),
                Property::string("schema").with_description("Schema (default: PUBLIC)"),
                Property::string("role").with_description("Role (optional)"),
            ],
            credentials_template: vec![
                Property::string("login").with_description("Username"),
                Property::string("password").as_password().with_description("Password"),
            ],
            aggregations: vec![],
            types_map,
        }
    }

    fn connection_params(&self, conn: &DataConnection) -> Result<ConnectionParams> {
        let p = Self::extract_params(conn)?;
        let mut options = HashMap::new();
        options.insert(OPT_ACCOUNT.to_string(), p.account);
        options.insert(OPT_USERNAME.to_string(), p.login);
        options.insert(OPT_PASSWORD.to_string(), p.password);
        if !p.db.is_empty() {
            options.insert(OPT_DB.to_string(), p.db);
        }
        if !p.warehouse.is_empty() {
            options.insert(OPT_WAREHOUSE.to_string(), p.warehouse);
        }
        if !p.schema.is_empty() {
            options.insert(OPT_SCHEMA.to_string(), p.schema);
        }
        if !p.role.is_empty() {
            options.insert(OPT_ROLE.to_string(), p.role);
        }
        options.insert(OPT_USE_HIGH_PRECISION.to_string(), "false".to_string());
        Ok(ConnectionParams { options })
    }

    fn test_connection(&self, conn: &DataConnection, pool_mgr: &ConnectionPoolManager) -> Result<()> {
        let params = Self::extract_params(conn)?;
        info!(account = %params.account, "Testing Snowflake connection");

        let conn_params = self.connection_params(conn)?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, conn, &conn_params)?;

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query("SELECT 1")?;
        let batches = stmt.execute_to_batches()?;
        if batches.is_empty() || batches[0].num_rows() < 1 {
            bail!("SELECT 1 returned no rows");
        }
        info!("Snowflake connection test successful");
        Ok(())
    }

    fn execute_query_streaming(
        &self,
        call: &FuncCall,
        pool_mgr: &ConnectionPoolManager,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
    ) -> Result<()> {
        let prepared = prepare_query(call, self)?;
        info!(sql = %prepared.sql, "Executing Snowflake query (streaming)");

        let conn_params = self.connection_params(call.connection())?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, call.connection(), &conn_params)?;

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query(&prepared.sql)?;
        stmt.execute_streaming(tx)?;
        info!("Snowflake streaming query completed");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::params::prepare_query;
    use crate::query::params::test_support::make_pattern_call;

    #[test]
    fn regex_pattern_uses_regexp_operator() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @email(email)",
            "email",
            "string",
            serde_json::json!({"op": "regex", "values": ["^\\w+@google.com.au$"]}),
        );
        let prepared = prepare_query(&call, &SnowflakeProvider).unwrap();
        assert_eq!(
            prepared.sql,
            "SELECT * FROM mock_data WHERE email REGEXP '^\\w+@google.com.au$'"
        );
    }
}
