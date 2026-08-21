use serde::{Deserialize, Serialize};

/// DataSource descriptor returned by GET /conn.
/// Mirrors Java: grok_connect/connectors_info/DataSource.java
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DataSourceDescriptor {
    #[serde(rename = "#type")]
    pub type_tag: String,
    #[serde(rename = "type")]
    pub source_type: String,
    /// Connector engine serving this type. Always `ADBC` here; the Java grok_connect
    /// omits it, so it stays optional on the wire.
    #[serde(default = "default_engine", skip_serializing_if = "Option::is_none")]
    pub engine: Option<String>,
    pub category: String,
    #[serde(default)]
    pub description: String,
    #[serde(default = "default_comment_start")]
    pub comment_start: String,
    #[serde(default = "default_query_language")]
    pub query_language: String,
    pub can_browse_schema: bool,
    #[serde(default)]
    pub default_schema: Option<String>,
    #[serde(default = "default_name_brackets")]
    pub name_brackets: String,
    #[serde(default = "default_true")]
    pub limit_at_end: bool,
    #[serde(default)]
    pub connection_template: Vec<Property>,
    #[serde(default)]
    pub credentials_template: Vec<Property>,
    #[serde(default)]
    pub aggregations: Vec<serde_json::Value>,
    #[serde(default)]
    pub types_map: std::collections::HashMap<String, String>,
}

pub const ENGINE_ADBC: &str = "ADBC";

pub fn default_engine() -> Option<String> { Some(ENGINE_ADBC.to_string()) }

/// The provider type without its engine suffix: `ClickHouseADBC` → `ClickHouse`.
///
/// The engine is carried by the descriptor's `engine` field, so the type name doesn't have
/// to spell it out as well — see the `stripEngineSuffix` option on `GET /conn` and `/info`.
pub fn strip_engine_suffix(source_type: &str) -> &str {
    match source_type.strip_suffix(ENGINE_ADBC) {
        Some(base) if !base.is_empty() => base,
        _ => source_type,
    }
}

fn default_comment_start() -> String { "--".to_string() }
fn default_query_language() -> String { "sql".to_string() }
fn default_name_brackets() -> String { "[]".to_string() }
fn default_true() -> bool { true }

/// Property definition for connection/credentials templates.
/// Mirrors Java: grok_connect/utils/Property.java
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct Property {
    pub property_type: String,
    pub name: String,
    // Omit empty optionals on serialize — Datlas' Property deserializer chokes on an
    // explicit `null` for these (its `$info` setter dereferences the value), and the Java
    // descriptors omit them too. `skip_serializing_if` keeps GET /conn parseable by Datlas.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub friendly_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub category: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub choices: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub info: Option<PropInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PropInfo {
    #[serde(default)]
    pub nullable: bool,
    #[serde(default)]
    pub password: bool,
    #[serde(default)]
    pub ignore: bool,
}

impl Property {
    pub fn string(name: &str) -> Self {
        Property {
            property_type: "string".to_string(),
            name: name.to_string(),
            friendly_name: None,
            description: None,
            category: None,
            format: None,
            choices: None,
            info: None,
        }
    }

    pub fn int(name: &str) -> Self {
        Property {
            property_type: "int".to_string(),
            name: name.to_string(),
            friendly_name: None,
            description: None,
            category: None,
            format: None,
            choices: None,
            info: None,
        }
    }

    pub fn bool(name: &str) -> Self {
        Property {
            property_type: "bool".to_string(),
            name: name.to_string(),
            friendly_name: None,
            description: None,
            category: None,
            format: None,
            choices: None,
            info: None,
        }
    }

    pub fn with_category(mut self, cat: &str) -> Self {
        self.category = Some(cat.to_string());
        self
    }

    pub fn with_description(mut self, desc: &str) -> Self {
        self.description = Some(desc.to_string());
        self
    }

    pub fn as_password(mut self) -> Self {
        self.info = Some(PropInfo { nullable: false, password: true, ignore: false });
        self
    }

    pub fn as_nullable(mut self) -> Self {
        self.info = Some(PropInfo { nullable: true, password: false, ignore: false });
        self
    }
}
