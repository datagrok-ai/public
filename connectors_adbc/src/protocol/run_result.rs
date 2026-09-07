use serde::{Deserialize, Serialize};

/// Query execution result/error envelope.
/// Mirrors Java: grok_connect/connectors_info/DataQueryRunResult.java
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct DataQueryRunResult {
    #[serde(rename = "#type")]
    pub type_tag: String,
    #[serde(default)]
    pub time_stamp: Option<String>,
    #[serde(default)]
    pub exec_time: f64,
    #[serde(default)]
    pub columns: i32,
    #[serde(default)]
    pub rows: i32,
    #[serde(default)]
    pub blob_length: i64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error_message: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error_stack_trace: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub log: Option<String>,
}

impl DataQueryRunResult {
    pub fn error(message: impl Into<String>) -> Self {
        DataQueryRunResult {
            type_tag: "DataQueryRunResult".to_string(),
            time_stamp: Some(chrono::Utc::now().to_rfc3339()),
            exec_time: 0.0,
            columns: 0,
            rows: 0,
            blob_length: 0,
            error_message: Some(message.into()),
            error_stack_trace: None,
            log: None,
        }
    }

    pub fn error_with_trace(message: impl Into<String>, trace: impl Into<String>) -> Self {
        let mut result = Self::error(message);
        result.error_stack_trace = Some(trace.into());
        result
    }

    pub fn success(columns: i32, rows: i32, exec_time: f64) -> Self {
        DataQueryRunResult {
            type_tag: "DataQueryRunResult".to_string(),
            time_stamp: Some(chrono::Utc::now().to_rfc3339()),
            exec_time,
            columns,
            rows,
            blob_length: 0,
            error_message: None,
            error_stack_trace: None,
            log: None,
        }
    }
}
