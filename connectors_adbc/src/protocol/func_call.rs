use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Top-level function call wrapper sent over WebSocket as "QUERY <json>".
/// Mirrors Java: grok_connect/connectors_info/FuncCall.java
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FuncCall {
    #[serde(default)]
    pub id: String,
    pub func: DataQuery,
    #[serde(default)]
    pub options: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub parameter_values: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub aux: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub log: Option<String>,
}

impl FuncCall {
    pub fn for_raw_sql(query: String, connection: DataConnection) -> FuncCall {
        serde_json::from_value(serde_json::json!({"func": {"query": query, "connection": connection}}))
            .expect("static FuncCall shape")
    }

    pub fn is_debug(&self) -> bool {
        self.options
            .get("debug")
            .and_then(|v| v.as_bool())
            .unwrap_or(false)
    }

    /// Get the DataConnection from the nested func.
    pub fn connection(&self) -> &DataConnection {
        &self.func.connection
    }

    /// Get the SQL query text.
    pub fn query(&self) -> &str {
        &self.func.query
    }
}

/// Query definition with SQL text, connection info, and parameters.
/// Mirrors Java: grok_connect/connectors_info/DataQuery.java
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DataQuery {
    #[serde(rename = "#type", default = "default_data_query_type")]
    pub type_tag: String,
    #[serde(default)]
    pub id: String,
    #[serde(default)]
    pub name: String,
    #[serde(default)]
    pub query: String,
    #[serde(default)]
    pub connection_id: String,
    pub connection: DataConnection,
    #[serde(default)]
    pub params: Vec<FuncParam>,
    #[serde(default)]
    pub options: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub aux: HashMap<String, serde_json::Value>,
}

fn default_data_query_type() -> String {
    "DataQuery".to_string()
}

/// Connection parameters and credentials.
/// Mirrors Java: grok_connect/connectors_info/DataConnection.java
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DataConnection {
    #[serde(default)]
    pub id: String,
    #[serde(default)]
    pub data_source: String,
    #[serde(default)]
    pub connection_string: Option<String>,
    #[serde(default)]
    pub credentials: Option<Credentials>,
    #[serde(default)]
    pub parameters: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub tags: HashMap<String, String>,
}

impl DataConnection {
    /// Get a string parameter (case-insensitive lookup).
    pub fn get_str(&self, key: &str) -> Option<&str> {
        let key_lower = key.to_lowercase();
        for (k, v) in &self.parameters {
            if k.to_lowercase() == key_lower {
                return v.as_str();
            }
        }
        None
    }

    /// Get an integer parameter.
    pub fn get_i64(&self, key: &str) -> Option<i64> {
        let key_lower = key.to_lowercase();
        for (k, v) in &self.parameters {
            if k.to_lowercase() == key_lower {
                if let Some(n) = v.as_i64() {
                    return Some(n);
                }
                if let Some(s) = v.as_str() {
                    return s.parse().ok();
                }
            }
        }
        None
    }

    /// Get credential parameter (login, password, token, etc.).
    pub fn get_credential(&self, key: &str) -> Option<&str> {
        self.credentials
            .as_ref()
            .and_then(|c| c.parameters.get(key))
            .and_then(|v| v.as_str())
    }

    /// Identity string for connection pooling keys — excludes secrets.
    pub fn connection_identity(&self) -> String {
        let host = self.get_str("server").or_else(|| self.get_str("host")).unwrap_or("");
        let port = self.get_i64("port").unwrap_or(0);
        let db = self.get_str("db").unwrap_or("");
        let path = self.get_str("httpPath").unwrap_or("");
        let warehouse = self.get_str("warehouse").unwrap_or("");
        let user = self.get_credential("login").unwrap_or("");
        format!("{}:{}:{}:{}:{}:{}", host, port, db, path, warehouse, user)
    }
}

/// Authentication credentials.
/// Mirrors Java: grok_connect/connectors_info/Credentials.java
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Credentials {
    #[serde(default)]
    pub parameters: HashMap<String, serde_json::Value>,
}

/// Function parameter definition.
/// Mirrors Java: grok_connect/connectors_info/FuncParam.java
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct FuncParam {
    #[serde(rename = "#type", default = "default_func_param_type")]
    pub type_tag: String,
    #[serde(default)]
    pub property_type: String,
    #[serde(default)]
    pub property_sub_type: Option<String>,
    #[serde(default)]
    pub name: String,
    #[serde(default)]
    pub value: serde_json::Value,
    #[serde(default)]
    pub is_input: bool,
    #[serde(default)]
    pub options: HashMap<String, String>,
}

fn default_func_param_type() -> String {
    "FuncParam".to_string()
}

/// Well-known parameter keys from DbCredentials.java.
#[allow(dead_code)]
pub mod db_credentials {
    pub const SERVER: &str = "server";
    pub const DB: &str = "db";
    pub const PORT: &str = "port";
    pub const LOGIN: &str = "login";
    pub const PASSWORD: &str = "password";
    pub const CONNECTION_STRING: &str = "connString";
    pub const SSL: &str = "ssl";
    pub const ACCESS_KEY: &str = "accessKey";
    pub const SECRET_KEY: &str = "secretKey";
}
