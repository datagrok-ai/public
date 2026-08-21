use crate::protocol::{default_engine, DataConnection, DataSourceDescriptor, FuncCall, Property};
use crate::providers::{Advertised, ConnectionParams, DbProvider};
use crate::query::params::{SqlDialect, prepare_query};
use crate::state::{AdbcConnection, ConnectionPoolManager, PoolDatabase};
use adbc_core::options::{AdbcVersion, OptionDatabase, OptionValue};
use adbc_core::{Driver as DriverTrait, Optionable};
use adbc_driver_manager::ManagedDriver;
use arrow::array::{Array, ArrayRef, StringArray, TimestampSecondArray, UInt16Array, UInt32Array};
use arrow::datatypes::{DataType, Field, Schema, TimeUnit};
use arrow::record_batch::RecordBatch;
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::sync::Arc;
use tracing::{info, warn};

// ClickHouse emits `LowCardinality(T)` as an Arrow `Dictionary` only when this session
// setting is on. It is NOT settable via ADBC `set_option` (driver whitelist) nor via a URL
// query param (the clickhouse-rs client strips it), so we issue it as a `SET` on the
// connection before each query (cheap; also re-arms it if a pooled session timed out).
const SET_DICT: &str = "SET output_format_arrow_low_cardinality_as_dictionary = 1";
// ClickHouse emits `String` as Arrow `Binary` unless this is on. The Datagrok Arrow decoder
// (`@datagrok-libraries/arrow` `fromFeather`) has no `Binary` case — it falls back to
// `Column.fromStrings`, which calls `.trim()` on the raw `Uint8Array` values and throws. Emitting
// `Utf8` makes strings decode as proper string columns (and lets dict-encoding kick in).
const SET_STRING_AS_STRING: &str = "SET output_format_arrow_string_as_string = 1";

const SECONDS_PER_DAY: i64 = 86_400;

pub struct ClickhouseProvider;

impl ClickhouseProvider {
    // No logging here: constructing a provider does no work, and `build_registry` logs
    // the probe that actually loads the driver.
    pub fn new() -> Self {
        ClickhouseProvider
    }

    /// Loaded as a C-ABI DLL via the driver manager: the driver's arrow 58 / adbc_core 0.23
    /// stays sealed behind the frozen C Data Interface, so we read its results (dictionary
    /// arrays included) with our arrow 57 / adbc_core 0.22.
    fn load_driver() -> Result<ManagedDriver> {
        match std::env::var("ADBC_CLICKHOUSE_DRIVER") {
            Ok(path) => ManagedDriver::load_dynamic_from_filename(
                &path, Some(b"AdbcClickhouseInit"), AdbcVersion::V110,
            ),
            Err(_) => ManagedDriver::load_dynamic_from_name(
                "adbc_clickhouse", Some(b"AdbcClickhouseInit"), AdbcVersion::V110,
            ),
        }
        .context(
            "Failed to load the ClickHouse ADBC driver. Build it with \
             'cargo build --release --features ffi' in the adbc_clickhouse crate and \
             drop adbc_clickhouse.dll into lib/, or set ADBC_CLICKHOUSE_DRIVER to its path."
        )
    }

    fn extract_params(conn: &DataConnection) -> Result<ClickhouseConnParams> {
        let host = conn.get_str("server")
            .or_else(|| conn.get_str("host"))
            .context("ClickHouse: server is required")?
            .to_string();
        let port = conn.get_i64("port").unwrap_or(8123) as u16;
        let login = conn.get_credential("login").unwrap_or("default").to_string();
        let password = conn.get_credential("password").unwrap_or("").to_string();
        Ok(ClickhouseConnParams { host, port, login, password })
    }

    fn prepare_session(adbc_conn: &mut AdbcConnection, conn: &DataConnection) -> Result<()> {
        for setting in [SET_DICT, SET_STRING_AS_STRING] {
            let mut stmt = adbc_conn.new_statement()?;
            stmt.set_sql_query(setting)?;
            stmt.execute_update()?;
        }
        Self::use_database(adbc_conn, conn)
    }

    /// Issue `USE <db>` when a non-default database is configured on the connection.
    fn use_database(adbc_conn: &mut AdbcConnection, conn: &DataConnection) -> Result<()> {
        let db = conn.get_str("db").unwrap_or("");
        if db.is_empty() || db.eq_ignore_ascii_case("default") {
            return Ok(());
        }
        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query(&format!("USE `{}`", db.replace('`', "``")))?;
        stmt.execute_update()?;
        Ok(())
    }

    /// The declared ClickHouse `(name, type)` of every column the query returns, in order.
    /// `DESCRIBE` infers the result schema without executing the query.
    /// Best-effort: on any failure we return nothing, which disables both temporal repair
    /// and UUID rewriting rather than failing the query.
    fn describe_columns(adbc_conn: &mut AdbcConnection, sql: &str) -> Vec<(String, String)> {
        match Self::describe(adbc_conn, sql) {
            Ok(columns) => columns,
            Err(e) => {
                warn!(error = %e, "DESCRIBE failed; temporal/UUID column repair disabled");
                Vec::new()
            }
        }
    }

    fn describe(adbc_conn: &mut AdbcConnection, sql: &str) -> Result<Vec<(String, String)>> {
        let inner = sql.trim().trim_end_matches(';');
        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query(&format!("DESCRIBE (\n{}\n)", inner))?;

        let mut columns = Vec::new();
        for batch in stmt.execute_to_batches()? {
            // DESCRIBE returns (name, type, default_type, ...); we want name + type.
            let names = batch
                .column(0)
                .as_any()
                .downcast_ref::<StringArray>()
                .context("DESCRIBE: `name` column is not a string")?;
            let types = batch
                .column(1)
                .as_any()
                .downcast_ref::<StringArray>()
                .context("DESCRIBE: `type` column is not a string")?;
            columns.extend((0..types.len())
                .map(|row| (names.value(row).to_string(), types.value(row).to_string())));
        }
        Ok(columns)
    }
}

/// The Datagrok Arrow decoder reads this field-metadata key; a value of [`BIGINT_METADATA_VALUE`]
/// tells it to rebuild the (String) column as an arbitrary-precision `bigint` instead of a
/// plain string. See [`annotate_bigint_columns`].
const BIGINT_METADATA_KEY: &str = "datagrok.type";
const BIGINT_METADATA_VALUE: &str = "bigint";

/// ClickHouse refuses to serialize some types into Arrow at all (`UUID`, `Code: 50`, even on
/// 26.x), and others arrive as a raw, sign-less `FixedSizeBinary` that `fromFeather` cannot
/// decode (`Int128`/`UInt128` -> 16 bytes, `Int256`/`UInt256` -> 32 bytes — the driver picks the
/// same Arrow type for the signed and unsigned pair, and also reuses it for IPv6, so there is no
/// way to recover the intended representation downstream of the query). Rewrite the query so
/// these columns come back as their canonical String form instead. `SELECT * REPLACE (expr AS
/// col)` swaps only those columns inside the `*` expansion, preserving column order and names, so
/// the declared-type indices used by [`repair_temporal_columns`] / [`annotate_bigint_columns`]
/// stay aligned. UUID matches the Java JDBC path, which also hands UUID back as a string. The
/// wide integers are additionally tagged as `bigint` (see [`annotate_bigint_columns`]) so they
/// decode to the same `bigint` column Java maps them to — no separate baseline needed.
/// No-op when the query returns none of these types (or when `DESCRIBE` was unavailable).
fn rewrite_unsupported_columns(sql: &str, columns: &[(String, String)]) -> String {
    let replacements: Vec<String> = columns
        .iter()
        .filter(|(_, declared)| {
            let declared = unwrap_nullable(declared);
            declared == "UUID" || is_wide_integer(declared)
        })
        .map(|(name, _)| {
            let ident = name.replace('`', "``");
            format!("toString(`{0}`) AS `{0}`", ident)
        })
        .collect();
    if replacements.is_empty() {
        return sql.to_string();
    }
    let inner = sql.trim().trim_end_matches(';');
    format!("SELECT * REPLACE ({}) FROM (\n{}\n)", replacements.join(", "), inner)
}

/// `Int128`/`Int256`/`UInt128`/`UInt256` — the wide integers that [`rewrite_unsupported_columns`]
/// casts to a decimal String, and that Java's JDBC connector maps to its `bigint` column type.
fn is_wide_integer(declared_type: &str) -> bool {
    matches!(declared_type, "Int128" | "Int256" | "UInt128" | "UInt256")
}

/// Tag the String columns produced from wide integers with `datagrok.type=bigint` field metadata,
/// so the Datagrok Arrow decoder rebuilds them as an arbitrary-precision `bigint` column (matching
/// Java) rather than a plain string. Column order/indices match the declared types because the
/// rewrite preserves them; a no-op when the query returns no wide integers.
fn annotate_bigint_columns(batch: RecordBatch, declared: &[String]) -> Result<RecordBatch> {
    let targets: Vec<usize> = (0..batch.num_columns())
        .filter(|index| {
            declared.get(*index).map(|t| is_wide_integer(unwrap_nullable(t))).unwrap_or(false)
        })
        .collect();
    if targets.is_empty() {
        return Ok(batch);
    }
    let mut fields: Vec<Field> =
        batch.schema().fields().iter().map(|f| f.as_ref().clone()).collect();
    for index in targets {
        let mut metadata = fields[index].metadata().clone();
        metadata.insert(BIGINT_METADATA_KEY.to_string(), BIGINT_METADATA_VALUE.to_string());
        fields[index] = fields[index].clone().with_metadata(metadata);
    }
    Ok(RecordBatch::try_new(Arc::new(Schema::new(fields)), batch.columns().to_vec())?)
}

/// ClickHouse's Arrow output drops the temporal type of the *narrow* date types: `Date`
/// arrives as `Uint16` (days since the epoch) and `DateTime` as `Uint32` (seconds). The wide
/// `Date32`/`DateTime64` already arrive as proper Arrow `Date32`/`Timestamp`, and a newer
/// ClickHouse may well emit the narrow ones correctly too. So a column is repaired only when
/// its declared type is temporal *and* it actually came through as the integer form — leaving
/// a genuine `UInt16` column, or a correctly-typed date, untouched.
///
/// Both are widened to `Timestamp(Second)` rather than `Date32`, so they land on the one
/// `fromFeather` branch that is known to decode to a Datagrok `datetime` column.
fn plan_temporal_repairs(schema: &Schema, declared: &[String]) -> Vec<(usize, i64)> {
    let mut plan = Vec::new();
    for (index, field) in schema.fields().iter().enumerate() {
        let Some(declared_type) = declared.get(index).map(|t| unwrap_nullable(t)) else {
            continue;
        };
        let seconds_per_unit = match (declared_type, field.data_type()) {
            ("Date", DataType::UInt16) => SECONDS_PER_DAY,
            (t, DataType::UInt32) if is_narrow_datetime(t) => 1,
            _ => continue,
        };
        plan.push((index, seconds_per_unit));
    }
    plan
}

fn unwrap_nullable(declared_type: &str) -> &str {
    declared_type
        .strip_prefix("Nullable(")
        .and_then(|inner| inner.strip_suffix(')'))
        .unwrap_or(declared_type)
}

/// `DateTime` or `DateTime('UTC')` — but never `DateTime64(...)`, which Arrow already gets right.
fn is_narrow_datetime(declared_type: &str) -> bool {
    declared_type == "DateTime" || declared_type.starts_with("DateTime(")
}

fn repair_temporal_columns(batch: RecordBatch, declared: &[String]) -> Result<RecordBatch> {
    let plan = plan_temporal_repairs(batch.schema_ref(), declared);
    if plan.is_empty() {
        return Ok(batch);
    }

    let mut columns: Vec<ArrayRef> = batch.columns().to_vec();
    let mut fields: Vec<Field> =
        batch.schema().fields().iter().map(|f| f.as_ref().clone()).collect();

    for (index, seconds_per_unit) in plan {
        let column = &columns[index];
        let seconds: TimestampSecondArray = match column.data_type() {
            DataType::UInt16 => column
                .as_any()
                .downcast_ref::<UInt16Array>()
                .context("expected a UInt16 Date column")?
                .iter()
                .map(|v| v.map(|days| days as i64 * seconds_per_unit))
                .collect(),
            DataType::UInt32 => column
                .as_any()
                .downcast_ref::<UInt32Array>()
                .context("expected a UInt32 DateTime column")?
                .iter()
                .map(|v| v.map(|seconds| seconds as i64 * seconds_per_unit))
                .collect(),
            other => bail!("unexpected temporal source type {}", other),
        };
        fields[index] = Field::new(
            fields[index].name(),
            DataType::Timestamp(TimeUnit::Second, None),
            fields[index].is_nullable(),
        );
        columns[index] = Arc::new(seconds) as ArrayRef;
    }

    Ok(RecordBatch::try_new(Arc::new(Schema::new(fields)), columns)?)
}

fn escape_literal(value: &str) -> String {
    value.replace('\\', "\\\\").replace('\'', "''")
}

#[derive(Debug)]
struct ClickhouseConnParams {
    host: String,
    port: u16,
    login: String,
    password: String,
}

/// The one provider this service currently serves.
impl Advertised for ClickhouseProvider {
    const ADVERTISED: bool = true;
}

/// ClickHouse spells regex matching `match(col, 'expr')` and needs `toDateTime` around
/// both sides of a datetime comparison; everything else is the vendor-neutral default.
impl SqlDialect for ClickhouseProvider {
    fn regex_match(&self, column: &str, expression: &str) -> Result<String> {
        Ok(format!("match({}, '{}')", column, expression.replace('\'', "''")))
    }

    fn datetime_operand(&self, operand: &str) -> String {
        format!("toDateTime({})", operand)
    }
}

impl DbProvider for ClickhouseProvider {
    fn name(&self) -> &str {
        "ClickHouseADBC"
    }

    fn probe(&self) -> Result<()> {
        Self::load_driver().map(|_| ())
    }

    fn build_pool_database(&self, params: &ConnectionParams) -> Result<PoolDatabase> {
        let mut driver = Self::load_driver()?;
        let mut db = DriverTrait::new_database(&mut driver)
            .context("Failed to create ClickHouse ADBC database")?;
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
                .with_context(|| format!("Failed to set ClickHouse option: {}", key_name))?;
        }
        Ok(Box::new(db))
    }

    fn descriptor(&self) -> DataSourceDescriptor {
        let mut types_map = HashMap::new();
        types_map.insert("uint8".to_string(), "int".to_string());
        types_map.insert("uint16".to_string(), "int".to_string());
        // UInt32's range overflows a signed int32, so Java maps it to bigint — mirror that.
        types_map.insert("uint32".to_string(), "bigint".to_string());
        types_map.insert("uint64".to_string(), "bigint".to_string());
        types_map.insert("int8".to_string(), "int".to_string());
        types_map.insert("int16".to_string(), "int".to_string());
        types_map.insert("int32".to_string(), "int".to_string());
        types_map.insert("int64".to_string(), "bigint".to_string());
        types_map.insert("int128".to_string(), "bigint".to_string());
        types_map.insert("int256".to_string(), "bigint".to_string());
        types_map.insert("uint128".to_string(), "bigint".to_string());
        types_map.insert("uint256".to_string(), "bigint".to_string());
        types_map.insert("float32".to_string(), "double".to_string());
        types_map.insert("float64".to_string(), "double".to_string());
        types_map.insert("#decimal.*".to_string(), "double".to_string());
        types_map.insert("string".to_string(), "string".to_string());
        types_map.insert("#fixedstring.*".to_string(), "string".to_string());
        types_map.insert("uuid".to_string(), "string".to_string());
        types_map.insert("date".to_string(), "datetime".to_string());
        types_map.insert("date32".to_string(), "datetime".to_string());
        types_map.insert("#datetime.*".to_string(), "datetime".to_string());
        types_map.insert("bool".to_string(), "bool".to_string());
        types_map.insert("#array.*".to_string(), "list".to_string());

        DataSourceDescriptor {
            type_tag: "DataSource".to_string(),
            // Distinct type so it coexists with the Java JDBC `ClickHouse` provider
            // instead of shadowing it; Datlas routes this type to this endpoint.
            source_type: "ClickHouseADBC".to_string(),
            engine: default_engine(),
            category: "Database".to_string(),
            description: "ClickHouse via ADBC (Arrow IPC, source-side dictionary)".to_string(),
            comment_start: "--".to_string(),
            query_language: "sql".to_string(),
            can_browse_schema: true,
            default_schema: Some("default".to_string()),
            name_brackets: "`".to_string(),
            limit_at_end: true,
            connection_template: vec![
                Property::string("server").with_description("ClickHouse host"),
                Property::int("port").with_description("HTTP port (default: 8123)"),
                Property::string("db").with_description("Database (default: default)"),
            ],
            credentials_template: vec![
                Property::string("login").with_description("Username (default: default)"),
                Property::string("password").as_password().with_description("Password"),
            ],
            aggregations: vec![],
            types_map,
        }
    }

    fn connection_params(&self, conn: &DataConnection) -> Result<ConnectionParams> {
        let p = Self::extract_params(conn)?;
        let mut options = HashMap::new();
        options.insert("uri".to_string(), format!("http://{}:{}", p.host, p.port));
        options.insert("username".to_string(), p.login);
        options.insert("password".to_string(), p.password);
        Ok(ConnectionParams { options })
    }

    fn test_connection(&self, conn: &DataConnection, pool_mgr: &ConnectionPoolManager) -> Result<()> {
        let params = Self::extract_params(conn)?;
        info!(host = %params.host, port = params.port, "Testing ClickHouse connection");

        let conn_params = self.connection_params(conn)?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, conn, &conn_params)?;

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query("SELECT 1")?;
        let batches = stmt.execute_to_batches()?;
        if batches.is_empty() || batches[0].num_rows() < 1 {
            bail!("SELECT 1 returned no rows");
        }
        info!("ClickHouse connection test successful");
        Ok(())
    }

    fn schemas_sql(&self) -> Option<String> {
        Some("SELECT DISTINCT table_schema FROM information_schema.columns".to_string())
    }

    fn schema_sql(&self, schema: Option<&str>, table: Option<&str>, include_key_info: bool) -> Option<String> {
        let mut sql = String::from(
            "SELECT c.table_schema AS table_schema, c.table_name AS table_name, \
             c.column_name AS column_name, c.data_type AS data_type, \
             if(t.table_type = 'VIEW', toInt8(1), toInt8(0)) AS is_view",
        );
        if include_key_info {
            sql.push_str(
                ", if(sc.is_in_primary_key = 1, toInt8(1), toInt8(0)) AS is_primary_key, \
                 toInt8(0) AS is_unique",
            );
        }
        sql.push_str(
            " FROM information_schema.columns c \
             JOIN information_schema.tables t \
               ON t.table_name = c.table_name \
              AND t.table_schema = c.table_schema \
              AND t.table_catalog = c.table_catalog",
        );
        if include_key_info {
            sql.push_str(
                " LEFT JOIN system.columns sc ON sc.database = c.table_catalog \
                 AND sc.table = c.table_name AND sc.name = c.column_name",
            );
        }
        if schema.is_some() || table.is_some() {
            sql.push_str(" WHERE");
            if let Some(s) = schema {
                sql.push_str(&format!(" c.table_schema = '{}'", escape_literal(s)));
            }
            if let Some(t) = table {
                if schema.is_some() {
                    sql.push_str(" AND");
                }
                sql.push_str(&format!(" c.table_name = '{}'", escape_literal(t)));
            }
        }
        Some(sql)
    }

    fn execute_query_streaming(
        &self,
        call: &FuncCall,
        pool_mgr: &ConnectionPoolManager,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
    ) -> Result<()> {
        let prepared = prepare_query(call, self)?;

        let conn_params = self.connection_params(call.connection())?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, call.connection(), &conn_params)?;
        Self::prepare_session(&mut adbc_conn, call.connection())?;
        let columns = Self::describe_columns(&mut adbc_conn, &prepared.sql);
        let declared: Vec<String> = columns.iter().map(|(_, t)| t.clone()).collect();
        let sql = rewrite_unsupported_columns(&prepared.sql, &columns);
        info!(sql = %sql, "Executing ClickHouse query (streaming)");

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query(&sql)?;
        stmt.execute_streaming_with(tx, |batch| repair_temporal_columns(batch, &declared)
            .and_then(|batch| annotate_bigint_columns(batch, &declared)))?;
        info!("ClickHouse streaming query completed");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query::params::prepare_query;
    use crate::query::params::test_support::make_pattern_call;

    fn schema(types: Vec<DataType>) -> Schema {
        Schema::new(
            types.into_iter()
                .enumerate()
                .map(|(i, t)| Field::new(format!("c{}", i), t, true))
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn narrow_date_and_datetime_are_repaired() {
        let schema = schema(vec![DataType::UInt16, DataType::UInt32]);
        let declared = vec!["Date".to_string(), "DateTime".to_string()];
        assert_eq!(plan_temporal_repairs(&schema, &declared), vec![(0, SECONDS_PER_DAY), (1, 1)]);
    }

    #[test]
    fn datetime_with_timezone_is_repaired() {
        let schema = schema(vec![DataType::UInt32]);
        let declared = vec!["DateTime('UTC')".to_string()];
        assert_eq!(plan_temporal_repairs(&schema, &declared), vec![(0, 1)]);
    }

    #[test]
    fn nullable_wrapper_is_unwrapped() {
        let schema = schema(vec![DataType::UInt16]);
        let declared = vec!["Nullable(Date)".to_string()];
        assert_eq!(plan_temporal_repairs(&schema, &declared), vec![(0, SECONDS_PER_DAY)]);
    }

    /// A newer ClickHouse already emits proper Arrow temporal types — leave them alone.
    #[test]
    fn already_correct_arrow_types_are_left_alone() {
        let schema = schema(vec![
            DataType::Date32,
            DataType::Timestamp(TimeUnit::Millisecond, None),
            DataType::Timestamp(TimeUnit::Second, None),
        ]);
        let declared = vec![
            "Date".to_string(),
            "DateTime64(3)".to_string(),
            "DateTime".to_string(),
        ];
        assert!(plan_temporal_repairs(&schema, &declared).is_empty());
    }

    /// A genuine integer column must never be turned into a timestamp.
    #[test]
    fn genuine_integer_columns_are_left_alone() {
        let schema = schema(vec![DataType::UInt16, DataType::UInt32, DataType::UInt64]);
        let declared = vec!["UInt16".to_string(), "UInt32".to_string(), "UInt64".to_string()];
        assert!(plan_temporal_repairs(&schema, &declared).is_empty());
    }

    /// DateTime64 is already a Timestamp; it must not be matched by the DateTime rule.
    #[test]
    fn datetime64_is_not_treated_as_narrow() {
        assert!(!is_narrow_datetime("DateTime64(3, 'Europe/Kyiv')"));
        assert!(is_narrow_datetime("DateTime"));
        assert!(is_narrow_datetime("DateTime('UTC')"));
    }

    /// Missing DESCRIBE output (best-effort fallback) must not repair anything.
    #[test]
    fn without_declared_types_nothing_is_repaired() {
        let schema = schema(vec![DataType::UInt16]);
        assert!(plan_temporal_repairs(&schema, &[]).is_empty());
    }

    #[test]
    fn uuid_columns_are_cast_to_string() {
        let cols = vec![
            ("id".to_string(), "UInt64".to_string()),
            ("uuid".to_string(), "UUID".to_string()),
            ("maybe".to_string(), "Nullable(UUID)".to_string()),
        ];
        assert_eq!(
            rewrite_unsupported_columns("SELECT * FROM uuid_type;", &cols),
            "SELECT * REPLACE (toString(`uuid`) AS `uuid`, toString(`maybe`) AS `maybe`) FROM (\nSELECT * FROM uuid_type\n)"
        );
    }

    #[test]
    fn wide_integer_columns_are_cast_to_string() {
        let cols = vec![
            ("int64_type".to_string(), "Int64".to_string()),
            ("int128_type".to_string(), "Int128".to_string()),
            ("int256_type".to_string(), "Int256".to_string()),
            ("uint128_type".to_string(), "UInt128".to_string()),
            ("uint256_type".to_string(), "UInt256".to_string()),
        ];
        assert_eq!(
            rewrite_unsupported_columns("SELECT * FROM t", &cols),
            "SELECT * REPLACE (toString(`int128_type`) AS `int128_type`, toString(`int256_type`) AS `int256_type`, toString(`uint128_type`) AS `uint128_type`, toString(`uint256_type`) AS `uint256_type`) FROM (\nSELECT * FROM t\n)"
        );
    }

    #[test]
    fn wide_integer_columns_are_tagged_bigint_but_uuid_is_not() {
        let schema = Schema::new(vec![
            Field::new("big", DataType::Utf8, true),
            Field::new("uuid", DataType::Utf8, true),
            Field::new("id", DataType::Int32, false),
        ]);
        let batch = RecordBatch::try_new(
            Arc::new(schema),
            vec![
                Arc::new(StringArray::from(vec!["1"])) as ArrayRef,
                Arc::new(StringArray::from(vec!["x"])) as ArrayRef,
                Arc::new(arrow::array::Int32Array::from(vec![1])) as ArrayRef,
            ],
        )
        .unwrap();
        let declared = vec!["UInt256".to_string(), "UUID".to_string(), "Int32".to_string()];
        let out = annotate_bigint_columns(batch, &declared).unwrap();
        assert_eq!(out.schema().field(0).metadata().get(BIGINT_METADATA_KEY), Some(&BIGINT_METADATA_VALUE.to_string()));
        assert!(out.schema().field(1).metadata().get(BIGINT_METADATA_KEY).is_none());
        assert!(out.schema().field(2).metadata().get(BIGINT_METADATA_KEY).is_none());
    }

    #[test]
    fn no_wide_integers_leaves_batch_untouched() {
        let schema = Schema::new(vec![Field::new("id", DataType::Int32, false)]);
        let batch = RecordBatch::try_new(
            Arc::new(schema),
            vec![Arc::new(arrow::array::Int32Array::from(vec![1])) as ArrayRef],
        )
        .unwrap();
        let out = annotate_bigint_columns(batch, &["Int32".to_string()]).unwrap();
        assert!(out.schema().field(0).metadata().is_empty());
    }

    #[test]
    fn queries_without_uuid_are_left_alone() {
        let cols = vec![("id".to_string(), "UInt64".to_string()), ("d".to_string(), "Date".to_string())];
        assert_eq!(rewrite_unsupported_columns("SELECT * FROM t", &cols), "SELECT * FROM t");
    }

    #[test]
    fn no_describe_output_means_no_rewrite() {
        assert_eq!(rewrite_unsupported_columns("SELECT * FROM t", &[]), "SELECT * FROM t");
    }

    #[test]
    fn regex_pattern_uses_match_function() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @email(email)",
            "email",
            "string",
            serde_json::json!({"op": "regex", "values": ["^\\w+@google.com.au$"]}),
        );
        let prepared = prepare_query(&call, &ClickhouseProvider).unwrap();
        assert_eq!(
            prepared.sql,
            "SELECT * FROM mock_data WHERE match(email, '^\\w+@google.com.au$')"
        );
    }

    #[test]
    fn datetime_pattern_wraps_both_operands() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @date(date)",
            "date",
            "datetime",
            serde_json::json!({"op": "before", "values": ["2022-01-01T00:00:00.000Z"], "include1": false}),
        );
        let prepared = prepare_query(&call, &ClickhouseProvider).unwrap();
        assert_eq!(
            prepared.sql,
            "SELECT * FROM mock_data WHERE (toDateTime(date) < toDateTime('2022-01-01 00:00:00'))"
        );
    }
}
