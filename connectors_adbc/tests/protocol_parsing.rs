use grok_connect_adbc::protocol::*;
use grok_connect_adbc::query::params::{DefaultDialect, prepare_query};
use grok_connect_adbc::providers::databricks::DatabricksProvider;

#[test]
fn test_parse_databricks_func_call() {
    let json = include_str!("../fixtures/func_call_databricks.json");
    let call: FuncCall = serde_json::from_str(json).unwrap();
    assert_eq!(call.id, "test_databricks_call");
    assert_eq!(call.func.type_tag, "DataQuery");
    assert_eq!(call.func.query, "SELECT * FROM samples.tpch.nation WHERE n_regionkey = @regionKey");
    assert_eq!(call.func.connection.data_source, "DatabricksADBC");
    assert_eq!(
        call.func.connection.get_str("server").unwrap(),
        "dbc-52240d8b-a70a.cloud.databricks.com"
    );
    assert_eq!(call.func.connection.get_i64("port").unwrap(), 443);
    assert_eq!(
        call.func.connection.get_str("httpPath").unwrap(),
        "/sql/1.0/warehouses/769f690b943edd7f"
    );
    assert_eq!(call.func.params.len(), 1);
    assert_eq!(call.func.params[0].name, "regionKey");
    assert_eq!(call.func.params[0].property_type, "int");
    assert_eq!(call.func.params[0].value, serde_json::json!(1));
}

#[test]
fn test_parse_snowflake_func_call() {
    let json = include_str!("../fixtures/func_call_snowflake.json");
    let call: FuncCall = serde_json::from_str(json).unwrap();
    assert_eq!(call.id, "test_snowflake_call");
    assert_eq!(call.func.connection.data_source, "SnowflakeADBC");
    assert_eq!(call.func.connection.get_str("server").unwrap(), "test_account.us-east-1");
    assert_eq!(call.func.connection.get_str("db").unwrap(), "TEST_DB");
    assert_eq!(call.func.connection.get_str("warehouse").unwrap(), "TEST_WH");
    assert_eq!(
        call.func.connection.get_credential("login").unwrap(),
        "test_user"
    );
    assert_eq!(call.func.params[0].name, "id");
}

#[test]
fn test_databricks_connection_parameters() {
    let json = include_str!("../fixtures/func_call_databricks.json");
    let call: FuncCall = serde_json::from_str(json).unwrap();
    let conn = call.connection();
    assert_eq!(conn.get_str("server").unwrap(), "dbc-52240d8b-a70a.cloud.databricks.com");
    assert_eq!(conn.get_i64("port").unwrap(), 443);
    assert_eq!(conn.get_str("httpPath").unwrap(), "/sql/1.0/warehouses/769f690b943edd7f");
    assert_eq!(conn.get_credential("password").unwrap(), "dapi_PLACEHOLDER_TOKEN");
}

#[test]
fn test_snowflake_connection_parameters() {
    let json = include_str!("../fixtures/func_call_snowflake.json");
    let call: FuncCall = serde_json::from_str(json).unwrap();
    let conn = call.connection();
    assert_eq!(conn.get_str("server").unwrap(), "test_account.us-east-1");
    assert_eq!(conn.get_str("db").unwrap(), "TEST_DB");
    assert_eq!(conn.get_str("warehouse").unwrap(), "TEST_WH");
    assert_eq!(conn.get_credential("login").unwrap(), "test_user");
    assert_eq!(conn.get_credential("password").unwrap(), "test_password");
}

#[test]
fn test_param_substitution_from_fixture() {
    let json = include_str!("../fixtures/func_call_databricks.json");
    let call: FuncCall = serde_json::from_str(json).unwrap();
    let prepared = prepare_query(&call, &DefaultDialect).unwrap();
    assert_eq!(
        prepared.sql,
        "SELECT * FROM samples.tpch.nation WHERE n_regionkey = 1"
    );
}

#[test]
fn test_jdbc_url_parsing() {
    let url = "jdbc:databricks://dbc-52240d8b-a70a.cloud.databricks.com:443/default;transportMode=http;ssl=1;AuthMech=3;httpPath=/sql/1.0/warehouses/769f690b943edd7f";
    let parts = DatabricksProvider::parse_jdbc_url(url).unwrap();
    assert_eq!(parts.host, "dbc-52240d8b-a70a.cloud.databricks.com");
    assert_eq!(parts.port, 443);
    assert_eq!(parts.http_path, "/sql/1.0/warehouses/769f690b943edd7f");
    assert_eq!(parts.catalog, "default");
}

#[test]
fn test_malformed_json_gives_error() {
    let bad_json = r#"{"not_valid_func_call": true}"#;
    let result = serde_json::from_str::<FuncCall>(bad_json);
    assert!(result.is_err());
}

#[test]
fn test_empty_query_connection() {
    let json = r##"{
        "id": "empty",
        "func": {
            "#type": "DataQuery",
            "query": "",
            "connection": {
                "dataSource": "Unknown",
                "parameters": {}
            },
            "params": []
        },
        "options": {},
        "parameterValues": {}
    }"##;
    let call: FuncCall = serde_json::from_str(json).unwrap();
    assert_eq!(call.func.query, "");
    assert_eq!(call.func.connection.data_source, "Unknown");
}

#[test]
fn test_connection_identity() {
    let json = include_str!("../fixtures/func_call_databricks.json");
    let call: FuncCall = serde_json::from_str(json).unwrap();
    let identity = call.connection().connection_identity();
    assert!(identity.contains("dbc-52240d8b-a70a.cloud.databricks.com"));
    assert!(identity.contains("443"));
}

#[test]
fn test_data_query_run_result_error() {
    let result = DataQueryRunResult::error("test error");
    assert_eq!(result.type_tag, "DataQueryRunResult");
    assert_eq!(result.error_message.as_deref(), Some("test error"));
    assert!(result.time_stamp.is_some());

    // Verify it serializes correctly
    let json = serde_json::to_string(&result).unwrap();
    assert!(json.contains("\"#type\":\"DataQueryRunResult\""));
    assert!(json.contains("\"errorMessage\":\"test error\""));
}

#[test]
fn test_data_query_run_result_success() {
    let result = DataQueryRunResult::success(5, 100, 1.5);
    assert_eq!(result.columns, 5);
    assert_eq!(result.rows, 100);
    assert_eq!(result.exec_time, 1.5);
    assert!(result.error_message.is_none());
}
