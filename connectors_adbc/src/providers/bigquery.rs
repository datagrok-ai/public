use crate::protocol::{default_engine, DataConnection, DataSourceDescriptor, FuncCall, Property};
use crate::providers::{Advertised, ConnectionParams, DbProvider};
use crate::query::params::{SqlDialect, prepare_query};
use crate::state::{ConnectionPoolManager, PoolDatabase};
use adbc_core::options::{AdbcVersion, OptionDatabase, OptionValue};
use adbc_core::{Driver as DriverTrait, Optionable};
use adbc_driver_manager::ManagedDriver;
use arrow::record_batch::RecordBatch;
use anyhow::{Context, Result, bail};
use base64::Engine;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use tracing::info;

// ADBC option keys for the BigQuery driver.
// See: apache/arrow-adbc go/adbc/driver/bigquery/driver.go.
const OPT_PROJECT_ID: &str = "adbc.bigquery.sql.project_id";
const OPT_DATASET_ID: &str = "adbc.bigquery.sql.dataset_id";
const OPT_AUTH_TYPE: &str = "adbc.bigquery.sql.auth_type";
const OPT_AUTH_CREDENTIALS: &str = "adbc.bigquery.sql.auth_credentials";
const OPT_AUTH_CLIENT_ID: &str = "adbc.bigquery.sql.auth.client_id";
const OPT_AUTH_CLIENT_SECRET: &str = "adbc.bigquery.sql.auth.client_secret";
const OPT_AUTH_REFRESH_TOKEN: &str = "adbc.bigquery.sql.auth.refresh_token";

// The driver validates `auth_type` against exactly these values and rejects anything else.
const AUTH_APPLICATION_DEFAULT: &str = "adbc.bigquery.sql.auth_type.auth_bigquery";
const AUTH_JSON_CREDENTIAL_STRING: &str = "adbc.bigquery.sql.auth_type.json_credential_string";
const AUTH_USER_AUTHENTICATION: &str = "adbc.bigquery.sql.auth_type.user_authentication";

// Credential field names, and the categories the platform groups them into. It sets
// `#chosen-auth-method` to the selected category — same contract as Java's
// `BigQueryDataProvider`, whose field names these reuse (`DbCredentials`).
const CHOSEN_AUTH_METHOD: &str = "#chosen-auth-method";
const METHOD_SERVICE_ACCOUNT: &str = "Service Account";
const METHOD_USER_AUTHENTICATION: &str = "User Authentication";
const CRED_SECRET_KEY: &str = "secretKey";
const CRED_CLIENT_ID: &str = "OAuth2ClientID";
const CRED_CLIENT_SECRET: &str = "OAuth2Secret";
const CRED_REFRESH_TOKEN: &str = "refreshToken";

/// The exported C-ABI entrypoint of the Go driver. Passed explicitly so that pointing
/// `ADBC_BIGQUERY_DRIVER` at a renamed copy still resolves — the driver manager would
/// otherwise derive the symbol from the file name.
const ENTRYPOINT: &[u8] = b"AdbcDriverBigqueryInit";

pub struct BigQueryProvider;

impl BigQueryProvider {
    pub fn new() -> Self {
        BigQueryProvider
    }

    /// Loaded as a C-ABI library via the driver manager, like ClickHouse's: the BigQuery
    /// driver is a Go build published only as a PyPI wheel — there is no crates.io crate to
    /// link, and (unlike Snowflake) no Debian package either.
    fn load_driver() -> Result<ManagedDriver> {
        match std::env::var("ADBC_BIGQUERY_DRIVER") {
            Ok(path) => ManagedDriver::load_dynamic_from_filename(
                &path, Some(ENTRYPOINT), AdbcVersion::V110,
            ),
            Err(_) => ManagedDriver::load_dynamic_from_name(
                "adbc_driver_bigquery", Some(ENTRYPOINT), AdbcVersion::V110,
            ),
        }
        .context(
            "Failed to load the BigQuery ADBC driver. Fetch it with \
             'scripts/fetch_adbc_driver.{sh,ps1} bigquery' from the grok_connect_adbc crate \
             root, or set ADBC_BIGQUERY_DRIVER to the full path of the library."
        )
    }

    fn extract_params(conn: &DataConnection) -> Result<BigQueryConnParams> {
        let project = conn.get_str("projectId")
            .filter(|value| !value.is_empty())
            .context("BigQuery: projectId is required")?
            .to_string();
        let dataset = conn.get_str("dataset").unwrap_or("").to_string();
        Ok(BigQueryConnParams { project, dataset, auth: Self::extract_auth(conn)? })
    }

    fn extract_auth(conn: &DataConnection) -> Result<BigQueryAuth> {
        let credential = |name: &str| conn.get_credential(name).filter(|value| !value.is_empty());
        match conn.get_credential(CHOSEN_AUTH_METHOD).unwrap_or("") {
            METHOD_SERVICE_ACCOUNT => Ok(BigQueryAuth::ServiceAccount {
                credentials: service_account_key(
                    credential(CRED_SECRET_KEY).context("BigQuery: secretKey is required")?)?,
            }),
            METHOD_USER_AUTHENTICATION => Ok(BigQueryAuth::UserAuthentication {
                client_id: credential(CRED_CLIENT_ID)
                    .context("BigQuery: OAuth2ClientID is required")?.to_string(),
                client_secret: credential(CRED_CLIENT_SECRET)
                    .context("BigQuery: OAuth2Secret is required")?.to_string(),
                refresh_token: credential(CRED_REFRESH_TOKEN)
                    .context("BigQuery: refreshToken is required")?.to_string(),
            }),
            // No method chosen: use a service-account key if one was supplied, and otherwise
            // fall back to the host's ambient identity rather than failing — that is how a
            // connection with no credentials at all is meant to work.
            "" => match credential(CRED_SECRET_KEY) {
                Some(key) => Ok(BigQueryAuth::ServiceAccount { credentials: service_account_key(key)? }),
                None => Ok(BigQueryAuth::ApplicationDefault),
            },
            other => bail!("BigQuery: unsupported authentication method '{}'", other),
        }
    }
}

#[derive(Debug)]
struct BigQueryConnParams {
    project: String,
    dataset: String,
    auth: BigQueryAuth,
}

/// How the driver authenticates. One variant per `auth_type` value the Go driver accepts.
#[derive(Debug, PartialEq)]
enum BigQueryAuth {
    /// A service-account key, handed to the driver as a JSON string.
    ServiceAccount { credentials: String },
    /// Installed-app OAuth — the driver exchanges the refresh token for access tokens itself.
    /// There is no option for handing it a ready-made access token, so the platform's generic
    /// OAuth flow (Java's `#token`) has no equivalent here.
    UserAuthentication { client_id: String, client_secret: String, refresh_token: String },
    /// Application Default Credentials: the ambient identity of the host — workload identity,
    /// the GCE metadata server, or `GOOGLE_APPLICATION_CREDENTIALS`.
    ApplicationDefault,
}

/// The service-account key in the form the driver wants it: raw JSON. A key pasted into the
/// credential field arrives verbatim, while one uploaded as a `.json` file arrives base64-encoded
/// (what Java's `BigQueryDataProvider` decodes) — accept both.
fn service_account_key(key: &str) -> Result<String> {
    let key = key.trim();
    if key.starts_with('{') {
        return Ok(key.to_string());
    }
    let decoded = base64::engine::general_purpose::STANDARD.decode(key)
        .context("BigQuery: secretKey is neither JSON nor base64")?;
    String::from_utf8(decoded).context("BigQuery: decoded secretKey is not valid UTF-8")
}

/// A non-cryptographic digest of the credentials on a connection, used only to keep pools of
/// otherwise identical connections apart. Never logged, and never leaves the process.
fn credentials_fingerprint(conn: &DataConnection) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    for name in [CRED_SECRET_KEY, CRED_CLIENT_ID, CRED_CLIENT_SECRET, CRED_REFRESH_TOKEN] {
        conn.get_credential(name).unwrap_or("").hash(&mut hasher);
    }
    hasher.finish()
}

/// Served: `BigQueryADBC` is a distinct type from the Java `BigQuery`, so the two coexist and
/// the user picks one when creating the connection. Requires the BigQuery ADBC driver to be
/// present, or the startup probe fails.
impl Advertised for BigQueryProvider {
    const ADVERTISED: bool = true;
}

/// BigQuery spells regex matching `REGEXP_CONTAINS(col, 'expr')`; the rest is the
/// vendor-neutral default. String literals take backslash escapes, so the expression — which
/// is a regex, and full of them — has to be escaped rather than only quote-doubled.
impl SqlDialect for BigQueryProvider {
    fn regex_match(&self, column: &str, expression: &str) -> Result<String> {
        Ok(format!(
            "REGEXP_CONTAINS({}, '{}')",
            column,
            expression.replace('\\', "\\\\").replace('\'', "\\'")
        ))
    }
}

impl DbProvider for BigQueryProvider {
    fn name(&self) -> &str {
        "BigQueryADBC"
    }

    fn probe(&self) -> Result<()> {
        Self::load_driver().map(|_| ())
    }

    fn build_pool_database(&self, params: &ConnectionParams) -> Result<PoolDatabase> {
        let mut driver = Self::load_driver()?;
        let mut db = DriverTrait::new_database(&mut driver)
            .context("Failed to create BigQuery ADBC database")?;
        // Every key BigQuery sends is a driver-specific `adbc.bigquery.sql.*` one — it has no
        // uri/username/password. `OptionDatabase`'s `From<&str>` still maps those three onto
        // their canonical variants, so a key added later lands in the right place by itself.
        for (key_name, value) in &params.options {
            db.set_option(OptionDatabase::from(key_name.as_str()), OptionValue::String(value.clone()))
                .with_context(|| format!("Failed to set BigQuery option: {}", key_name))?;
        }
        Ok(Box::new(db))
    }

    /// `connection_identity` keys on host/port/db/user, none of which a BigQuery connection
    /// carries — every one of them would land in the same pool and reuse whichever credentials
    /// opened it first. Key on the project and dataset plus a fingerprint of the credentials.
    fn pool_identity(&self, conn: &DataConnection) -> String {
        format!(
            "{}:{}:{:016x}",
            conn.get_str("projectId").unwrap_or(""),
            conn.get_str("dataset").unwrap_or(""),
            credentials_fingerprint(conn),
        )
    }

    fn descriptor(&self) -> DataSourceDescriptor {
        let mut types_map = HashMap::new();
        // Mirrors Java BigQueryDataProvider.typesMap.
        types_map.insert("string".to_string(), "string".to_string());
        types_map.insert("#int.*".to_string(), "int".to_string());
        types_map.insert("#float.*".to_string(), "double".to_string());
        types_map.insert("numeric".to_string(), "double".to_string());
        types_map.insert("bignumeric".to_string(), "bigint".to_string());
        types_map.insert("boolean".to_string(), "bool".to_string());
        types_map.insert("date".to_string(), "datetime".to_string());
        types_map.insert("datetime".to_string(), "datetime".to_string());
        types_map.insert("time".to_string(), "datetime".to_string());
        types_map.insert("#timestamp.*".to_string(), "datetime".to_string());
        types_map.insert("#interval.*".to_string(), "string".to_string());
        types_map.insert("array".to_string(), "list".to_string());
        types_map.insert("bytes".to_string(), "blob".to_string());

        DataSourceDescriptor {
            type_tag: "DataSource".to_string(),
            // Distinct type so it coexists with the Java JDBC `BigQuery` provider instead of
            // shadowing it; Datlas routes this type to this endpoint.
            source_type: "BigQueryADBC".to_string(),
            engine: default_engine(),
            category: "Database".to_string(),
            description: "BigQuery via ADBC (Arrow IPC)".to_string(),
            comment_start: "--".to_string(),
            query_language: "sql".to_string(),
            can_browse_schema: false,
            default_schema: None,
            name_brackets: "`".to_string(),
            limit_at_end: true,
            connection_template: vec![
                Property::string("projectId").with_description("GCP project ID (billing project)"),
                Property::string("dataset")
                    .with_description("Default dataset for unqualified table names (optional)"),
            ],
            credentials_template: vec![
                Property::string(CRED_SECRET_KEY)
                    .as_password()
                    .with_category(METHOD_SERVICE_ACCOUNT)
                    .with_description("Service account key, as the .json file's contents or base64 of it"),
                Property::string(CRED_CLIENT_ID)
                    .with_category(METHOD_USER_AUTHENTICATION)
                    .with_description("OAuth client ID"),
                Property::string(CRED_CLIENT_SECRET)
                    .as_password()
                    .with_category(METHOD_USER_AUTHENTICATION)
                    .with_description("OAuth client secret"),
                Property::string(CRED_REFRESH_TOKEN)
                    .as_password()
                    .with_category(METHOD_USER_AUTHENTICATION)
                    .with_description("OAuth refresh token"),
            ],
            aggregations: vec![],
            types_map,
        }
    }

    fn connection_params(&self, conn: &DataConnection) -> Result<ConnectionParams> {
        let p = Self::extract_params(conn)?;
        let mut options = HashMap::new();
        options.insert(OPT_PROJECT_ID.to_string(), p.project);
        if !p.dataset.is_empty() {
            options.insert(OPT_DATASET_ID.to_string(), p.dataset);
        }
        match p.auth {
            BigQueryAuth::ServiceAccount { credentials } => {
                options.insert(OPT_AUTH_TYPE.to_string(), AUTH_JSON_CREDENTIAL_STRING.to_string());
                options.insert(OPT_AUTH_CREDENTIALS.to_string(), credentials);
            }
            BigQueryAuth::UserAuthentication { client_id, client_secret, refresh_token } => {
                options.insert(OPT_AUTH_TYPE.to_string(), AUTH_USER_AUTHENTICATION.to_string());
                options.insert(OPT_AUTH_CLIENT_ID.to_string(), client_id);
                options.insert(OPT_AUTH_CLIENT_SECRET.to_string(), client_secret);
                options.insert(OPT_AUTH_REFRESH_TOKEN.to_string(), refresh_token);
            }
            BigQueryAuth::ApplicationDefault => {
                options.insert(OPT_AUTH_TYPE.to_string(), AUTH_APPLICATION_DEFAULT.to_string());
            }
        }
        Ok(ConnectionParams { options })
    }

    fn test_connection(&self, conn: &DataConnection, pool_mgr: &ConnectionPoolManager) -> Result<()> {
        let params = Self::extract_params(conn)?;
        info!(project = %params.project, "Testing BigQuery connection");

        let conn_params = self.connection_params(conn)?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, conn, &conn_params)?;

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query("SELECT 1")?;
        let batches = stmt.execute_to_batches()?;
        if batches.is_empty() || batches[0].num_rows() < 1 {
            bail!("SELECT 1 returned no rows");
        }
        info!("BigQuery connection test successful");
        Ok(())
    }

    fn execute_query_streaming(
        &self,
        call: &FuncCall,
        pool_mgr: &ConnectionPoolManager,
        tx: tokio::sync::mpsc::Sender<RecordBatch>,
    ) -> Result<()> {
        let prepared = prepare_query(call, self)?;
        info!(sql = %prepared.sql, "Executing BigQuery query (streaming)");

        let conn_params = self.connection_params(call.connection())?;
        let mut adbc_conn = pool_mgr.get_or_create_connection(self, call.connection(), &conn_params)?;

        let mut stmt = adbc_conn.new_statement()?;
        stmt.set_sql_query(&prepared.sql)?;
        stmt.execute_streaming(tx)?;
        info!("BigQuery streaming query completed");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protocol::Credentials;
    use crate::query::params::test_support::make_pattern_call;

    fn connection(params: Vec<(&str, &str)>, credentials: Vec<(&str, &str)>) -> DataConnection {
        let to_map = |pairs: Vec<(&str, &str)>| pairs
            .into_iter()
            .map(|(k, v)| (k.to_string(), serde_json::Value::String(v.to_string())))
            .collect();
        DataConnection {
            id: String::new(),
            data_source: "BigQueryADBC".to_string(),
            connection_string: None,
            credentials: Some(Credentials { parameters: to_map(credentials) }),
            parameters: to_map(params),
            tags: HashMap::new(),
        }
    }

    #[test]
    fn regex_pattern_uses_regexp_contains() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @email(email)",
            "email",
            "string",
            serde_json::json!({"op": "regex", "values": ["^\\w+@google.com.au$"]}),
        );
        let prepared = prepare_query(&call, &BigQueryProvider).unwrap();
        assert_eq!(
            prepared.sql,
            "SELECT * FROM mock_data WHERE REGEXP_CONTAINS(email, '^\\\\w+@google.com.au$')"
        );
    }

    #[test]
    fn service_account_key_accepts_raw_json_and_base64() {
        let json = r#"{"type":"service_account"}"#;
        let encoded = base64::engine::general_purpose::STANDARD.encode(json);
        assert_eq!(service_account_key(json).unwrap(), json);
        assert_eq!(service_account_key(&encoded).unwrap(), json);
    }

    #[test]
    fn service_account_key_rejects_garbage() {
        assert!(service_account_key("not json, not base64!!").is_err());
    }

    #[test]
    fn chosen_auth_method_selects_the_credential_set() {
        let conn = connection(
            vec![("projectId", "p")],
            vec![
                (CHOSEN_AUTH_METHOD, METHOD_USER_AUTHENTICATION),
                (CRED_SECRET_KEY, r#"{"type":"service_account"}"#),
                (CRED_CLIENT_ID, "id"),
                (CRED_CLIENT_SECRET, "secret"),
                (CRED_REFRESH_TOKEN, "refresh"),
            ],
        );
        assert_eq!(
            BigQueryProvider::extract_auth(&conn).unwrap(),
            BigQueryAuth::UserAuthentication {
                client_id: "id".to_string(),
                client_secret: "secret".to_string(),
                refresh_token: "refresh".to_string(),
            }
        );
    }

    /// A connection carrying no credentials falls back to the host's ambient identity.
    #[test]
    fn no_credentials_means_application_default() {
        let conn = connection(vec![("projectId", "p")], vec![]);
        assert_eq!(BigQueryProvider::extract_auth(&conn).unwrap(), BigQueryAuth::ApplicationDefault);
    }

    #[test]
    fn service_account_key_alone_selects_service_account() {
        let json = r#"{"type":"service_account"}"#;
        let conn = connection(vec![("projectId", "p")], vec![(CRED_SECRET_KEY, json)]);
        assert_eq!(
            BigQueryProvider::extract_auth(&conn).unwrap(),
            BigQueryAuth::ServiceAccount { credentials: json.to_string() }
        );
    }

    #[test]
    fn missing_project_id_is_an_error() {
        let conn = connection(vec![], vec![]);
        assert!(BigQueryProvider::extract_params(&conn).is_err());
    }

    #[test]
    fn auth_type_and_project_reach_the_adbc_options() {
        let json = r#"{"type":"service_account"}"#;
        let conn = connection(
            vec![("projectId", "my-project"), ("dataset", "my_dataset")],
            vec![(CHOSEN_AUTH_METHOD, METHOD_SERVICE_ACCOUNT), (CRED_SECRET_KEY, json)],
        );
        let options = BigQueryProvider.connection_params(&conn).unwrap().options;
        assert_eq!(options[OPT_PROJECT_ID], "my-project");
        assert_eq!(options[OPT_DATASET_ID], "my_dataset");
        assert_eq!(options[OPT_AUTH_TYPE], AUTH_JSON_CREDENTIAL_STRING);
        assert_eq!(options[OPT_AUTH_CREDENTIALS], json);
    }

    /// Two connections to the same project under different keys must not share a pool —
    /// they would share ADBC connections, and so the credentials behind them.
    #[test]
    fn pool_identity_separates_different_credentials() {
        let project = vec![("projectId", "p")];
        let one = connection(project.clone(), vec![(CRED_SECRET_KEY, r#"{"a":1}"#)]);
        let two = connection(project.clone(), vec![(CRED_SECRET_KEY, r#"{"b":2}"#)]);
        assert_ne!(
            BigQueryProvider.pool_identity(&one),
            BigQueryProvider.pool_identity(&two)
        );
        assert_eq!(
            BigQueryProvider.pool_identity(&one),
            BigQueryProvider.pool_identity(&connection(project, vec![(CRED_SECRET_KEY, r#"{"a":1}"#)]))
        );
    }

    #[test]
    fn pool_identity_separates_different_projects() {
        let credentials = vec![(CRED_SECRET_KEY, r#"{"a":1}"#)];
        assert_ne!(
            BigQueryProvider.pool_identity(&connection(vec![("projectId", "one")], credentials.clone())),
            BigQueryProvider.pool_identity(&connection(vec![("projectId", "two")], credentials))
        );
    }
}
