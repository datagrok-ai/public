use crate::protocol::{FuncCall, FuncParam};
use anyhow::{Result, bail};
use regex::Regex;

/// Result of parameter substitution: the rewritten SQL and ordered parameter values.
#[derive(Debug)]
pub struct PreparedQuery {
    pub sql: String,
    pub values: Vec<serde_json::Value>,
}

/// Expand Datagrok search patterns (`@id(id)`) into SQL predicates and replace plain
/// `@paramName` placeholders with SQL literals.
///
/// Both happen in a **single** regex pass. `replace_all` never re-scans what it just
/// substituted, so an inlined value that itself contains `@` — an email regex, say —
/// cannot be mistaken for another parameter.
///
/// Fails when the dialect cannot render a pattern the user asked for; the query must not
/// run with that filter quietly dropped.
pub fn prepare_query(call: &FuncCall, dialect: &dyn SqlDialect) -> Result<PreparedQuery> {
    let sql = &call.func.query;
    let params = &call.func.params;
    let param_values = &call.parameter_values;
    let patterns = call.options.get("patterns");

    if params.is_empty() && param_values.is_empty() && patterns.is_none() {
        return Ok(PreparedQuery { sql: sql.clone(), values: Vec::new() });
    }

    let re = Regex::new(r"@(\w+)(\((\S+?)\))?").unwrap();
    let mut values = Vec::new();
    // `replace_all`'s closure cannot fail, so the first rendering error is parked here and
    // returned once the pass is over.
    let mut failure = None;

    let result = re.replace_all(sql, |caps: &regex::Captures| {
        let name = &caps[1];
        if let Some(col_name) = caps.get(3) {
            match pattern_to_sql(call, patterns, name, col_name.as_str(), dialect) {
                Ok(Some(predicate)) => return predicate,
                Ok(None) => {}
                Err(e) => {
                    failure.get_or_insert(e);
                    return String::new();
                }
            }
        }
        let value = resolve_param_value(name, params, param_values);
        let replacement = format_value_for_sql(&value);
        values.push(value);
        // Not a search pattern — keep whatever parenthesized text followed the parameter.
        match caps.get(2) {
            Some(suffix) => format!("{}{}", replacement, suffix.as_str()),
            None => replacement,
        }
    });

    match failure {
        Some(e) => Err(e),
        None => Ok(PreparedQuery { sql: result.into_owned(), values }),
    }
}

/// Resolve a parameter value: check parameterValues first, then params list.
fn resolve_param_value(
    name: &str,
    params: &[FuncParam],
    param_values: &std::collections::HashMap<String, serde_json::Value>,
) -> serde_json::Value {
    // Check parameterValues map first (case-insensitive)
    let name_lower = name.to_lowercase();
    for (k, v) in param_values {
        if k.to_lowercase() == name_lower {
            return v.clone();
        }
    }
    // Fall back to params list
    for p in params {
        if p.name.to_lowercase() == name_lower {
            return p.value.clone();
        }
    }
    serde_json::Value::Null
}

/// Format a JSON value as a SQL literal for inline substitution.
pub fn format_value_for_sql(value: &serde_json::Value) -> String {
    match value {
        serde_json::Value::Null => "NULL".to_string(),
        serde_json::Value::Bool(b) => if *b { "TRUE" } else { "FALSE" }.to_string(),
        serde_json::Value::Number(n) => n.to_string(),
        serde_json::Value::String(s) => {
            // Escape single quotes for SQL safety
            let escaped = s.replace('\'', "''");
            format!("'{}'", escaped)
        }
        _ => format!("'{}'", value),
    }
}

// Search-pattern operators, as emitted by the platform into `FuncCall.options["patterns"]`.
// Mirrors grok_connect's `utils/PatternMatcher.java`.
pub const OP_NONE: &str = "none";
pub const OP_RANGE_NUM: &str = "-";
pub const OP_RANGE_DATE_TIME: &str = "range";
pub const OP_CONTAINS: &str = "contains";
pub const OP_STARTS_WITH: &str = "starts with";
pub const OP_ENDS_WITH: &str = "ends with";
pub const OP_EQUALS: &str = "equals";
pub const OP_REGEXP: &str = "regex";
pub const OP_IN: &str = "in";
pub const OP_NOT_IN: &str = "not in";
pub const OP_BEFORE: &str = "before";
pub const OP_AFTER: &str = "after";
pub const OP_IS_NULL: &str = "is null";
pub const OP_IS_NOT_NULL: &str = "is not null";

pub const TRUE_PREDICATE: &str = "(1 = 1)";

/// A search-pattern matcher, already parsed by the platform.
pub struct PatternMatcher {
    pub op: String,
    pub values: Vec<serde_json::Value>,
    pub include1: bool,
    pub include2: bool,
}

impl PatternMatcher {
    fn from_json(value: &serde_json::Value) -> Option<Self> {
        Some(PatternMatcher {
            op: value.get("op")?.as_str()?.to_string(),
            values: value.get("values").and_then(|v| v.as_array()).cloned().unwrap_or_default(),
            include1: value.get("include1").and_then(|v| v.as_bool()).unwrap_or(false),
            include2: value.get("include2").and_then(|v| v.as_bool()).unwrap_or(false),
        })
    }

    pub fn value(&self, index: usize) -> Option<&serde_json::Value> {
        self.values.get(index)
    }

    pub fn raw_string_value(&self, index: usize) -> String {
        match self.value(index) {
            Some(serde_json::Value::String(s)) => s.clone(),
            Some(other) => other.to_string(),
            None => String::new(),
        }
    }

    /// `LIKE` comparisons are case-insensitive in grok_connect, which lowercases both sides.
    pub fn lowercase_value(&self, index: usize) -> String {
        self.raw_string_value(index).to_lowercase()
    }
}

/// Comparison operator for `before`/`after`, honouring the inclusive flag.
pub fn cmp(op: &str, include: bool) -> &'static str {
    match op {
        OP_BEFORE => if include { " <= " } else { " < " },
        OP_AFTER => if include { " >= " } else { " > " },
        _ => "",
    }
}

/// How a provider renders search patterns as SQL.
///
/// Every method has a default that emits vendor-neutral SQL, so a provider implements
/// this trait in **its own** `providers/<name>.rs` and overrides only what its dialect
/// spells differently — the Rust counterpart of the overridable hooks in Java's
/// `JdbcDataProvider`. Nothing here needs editing to add a provider.
///
/// Three levels of granularity, coarser ones built out of finer ones:
/// - operand hooks (`regex_match`, `datetime_operand`, `like`) — the usual case;
/// - per-type hooks (`numeric_pattern`, `string_pattern`, …) — replace a whole family;
/// - `pattern` — the dispatcher, for providers that add pattern types of their own.
pub trait SqlDialect {
    /// A regex-match predicate.
    ///
    /// The default errors, mirroring Java's `JdbcDataProvider.getRegexQuery`, which throws
    /// `UnsupportedOperationException`. A dialect that cannot express the filter the user
    /// asked for must say so — degrading to [`TRUE_PREDICATE`] would hand back the whole
    /// unfiltered table as if it were the filtered result.
    fn regex_match(&self, _column: &str, _expression: &str) -> Result<String> {
        bail!("REGEXP is not supported for this provider")
    }

    /// Wraps a datetime operand (column or literal) before comparison.
    fn datetime_operand(&self, operand: &str) -> String {
        operand.to_string()
    }

    /// A case-insensitive `LIKE`. `pattern` is already lowercased and wildcarded.
    fn like(&self, column: &str, pattern: &str) -> String {
        format!("(LOWER({}) LIKE {})", column, quote(pattern))
    }

    fn numeric_pattern(&self, column: &str, matcher: &PatternMatcher) -> Result<String> {
        Ok(match matcher.op.as_str() {
            OP_NONE => TRUE_PREDICATE.to_string(),
            OP_RANGE_NUM => format!(
                "({} >= {} AND {} <= {})",
                column,
                format_value_at(matcher, 0),
                column,
                format_value_at(matcher, 1)
            ),
            OP_IN | OP_NOT_IN => in_predicate(column, matcher),
            OP_IS_NULL | OP_IS_NOT_NULL => format!("({} {})", column, matcher.op),
            // Anything else is the comparison operator itself: `>`, `>=`, `<`, `<=`, `=`, `!=`.
            op => format!("({} {} {})", column, op, format_value_at(matcher, 0)),
        })
    }

    fn bigint_pattern(&self, column: &str, matcher: &PatternMatcher) -> Result<String> {
        Ok(match matcher.op.as_str() {
            OP_EQUALS => format!("({} = {})", column, format_value_at(matcher, 0)),
            OP_IN | OP_NOT_IN => in_predicate(column, matcher),
            _ => TRUE_PREDICATE.to_string(),
        })
    }

    fn string_pattern(&self, column: &str, matcher: &PatternMatcher) -> Result<String> {
        let value = matcher.lowercase_value(0);
        Ok(match matcher.op.as_str() {
            OP_EQUALS => self.like(column, &value),
            OP_CONTAINS => self.like(column, &format!("%{}%", value)),
            OP_STARTS_WITH => self.like(column, &format!("{}%", value)),
            OP_ENDS_WITH => self.like(column, &format!("%{}", value)),
            OP_REGEXP => self.regex_match(column, &value)?,
            OP_IN | OP_NOT_IN => in_predicate(column, matcher),
            OP_IS_NULL | OP_IS_NOT_NULL => format!("({} {})", column, matcher.op),
            _ => TRUE_PREDICATE.to_string(),
        })
    }

    fn datetime_pattern(&self, column: &str, matcher: &PatternMatcher) -> Result<String> {
        let operand = self.datetime_operand(column);
        let literal = |index: usize| {
            self.datetime_operand(&quote(&normalize_datetime(&matcher.raw_string_value(index))))
        };

        Ok(match matcher.op.as_str() {
            OP_EQUALS => format!("({} = {})", operand, literal(0)),
            OP_BEFORE | OP_AFTER => {
                format!("({}{}{})", operand, cmp(&matcher.op, matcher.include1), literal(0))
            }
            OP_RANGE_DATE_TIME => format!(
                "({}{}{} AND {}{}{})",
                operand,
                cmp(OP_AFTER, matcher.include1),
                literal(0),
                operand,
                cmp(OP_BEFORE, matcher.include2),
                literal(1)
            ),
            OP_IS_NULL | OP_IS_NOT_NULL => format!("({} {})", column, matcher.op),
            _ => TRUE_PREDICATE.to_string(),
        })
    }

    fn bool_pattern(&self, column: &str, matcher: &PatternMatcher) -> Result<String> {
        if matcher.op != OP_EQUALS {
            return Ok(TRUE_PREDICATE.to_string());
        }
        let value = matcher.value(0).cloned().unwrap_or(serde_json::Value::Bool(true));
        Ok(format!("({} = {})", column, format_value_for_sql(&value)))
    }

    /// Dispatch on the pattern type declared in the parameter's `pattern` option.
    fn pattern(&self, pattern_type: &str, column: &str, matcher: &PatternMatcher) -> Result<String> {
        match pattern_type {
            "num" | "double" | "int" => self.numeric_pattern(column, matcher),
            "bigint" => self.bigint_pattern(column, matcher),
            "string" => self.string_pattern(column, matcher),
            "datetime" => self.datetime_pattern(column, matcher),
            "bool" => self.bool_pattern(column, matcher),
            _ => Ok(TRUE_PREDICATE.to_string()),
        }
    }
}

/// Vendor-neutral rendering — the dialect of a provider that overrides nothing.
pub struct DefaultDialect;

impl SqlDialect for DefaultDialect {}

/// The SQL predicate a `@paramName(columnName)` search pattern denotes.
/// `Ok(None)` when the parameter is not a search pattern — the caller then treats it as a plain one.
fn pattern_to_sql(
    call: &FuncCall,
    patterns: Option<&serde_json::Value>,
    param_name: &str,
    col_name: &str,
    dialect: &dyn SqlDialect,
) -> Result<Option<String>> {
    let Some(param) = find_param(param_name, &call.func.params) else {
        return Ok(None);
    };
    // No `pattern` option: not a search pattern, leave it to the plain-parameter pass.
    let Some(pattern_type) = param.options.get("pattern") else {
        return Ok(None);
    };

    // A pattern parameter with no matcher applies no filter. Never fall through to the
    // plain-parameter pass here — that would inline the raw expression and emit `'>=29'(id)`.
    let matcher = patterns
        .and_then(|p| p.get(&param.name))
        .and_then(PatternMatcher::from_json);
    let Some(matcher) = matcher else {
        return Ok(Some(TRUE_PREDICATE.to_string()));
    };

    dialect.pattern(pattern_type, col_name, &matcher).map(Some)
}

fn find_param<'a>(name: &str, params: &'a [FuncParam]) -> Option<&'a FuncParam> {
    params.iter().find(|p| p.name.eq_ignore_ascii_case(name))
}

/// `(col in (1, 2))` / `(col not in ('a', 'b'))`
pub fn in_predicate(col_name: &str, matcher: &PatternMatcher) -> String {
    let values: Vec<String> = matcher.values.iter().map(format_value_for_sql).collect();
    format!("({} {} ({}))", col_name, matcher.op, values.join(", "))
}

pub fn format_value_at(matcher: &PatternMatcher, index: usize) -> String {
    match matcher.value(index) {
        Some(value) => format_value_for_sql(value),
        None => "NULL".to_string(),
    }
}

pub fn quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', "''"))
}

/// The platform sends ISO-8601 (`2022-01-01T00:00:00.000Z`); SQL wants `2022-01-01 00:00:00`.
pub fn normalize_datetime(value: &str) -> String {
    let value = value.replace('T', " ");
    let cut = value.find('.').or_else(|| value.find('+')).unwrap_or(value.len());
    value[..cut].trim_end_matches(|c| c == 'Z' || c == 'z').trim_end().to_string()
}

/// FuncCall builders shared with the per-provider dialect tests in `providers/*.rs`.
#[cfg(test)]
pub mod test_support {
    use crate::protocol::{DataConnection, DataQuery, FuncCall, FuncParam};
    use std::collections::HashMap;

    pub fn make_call(query: &str, params: Vec<(&str, &str, serde_json::Value)>) -> FuncCall {
        let func_params: Vec<FuncParam> = params
            .iter()
            .map(|(name, ptype, val)| FuncParam {
                type_tag: "FuncParam".to_string(),
                property_type: ptype.to_string(),
                property_sub_type: None,
                name: name.to_string(),
                value: val.clone(),
                is_input: true,
                options: HashMap::new(),
            })
            .collect();

        let param_values: HashMap<String, serde_json::Value> = params
            .into_iter()
            .map(|(name, _, val)| (name.to_string(), val))
            .collect();

        FuncCall {
            id: "test".to_string(),
            func: DataQuery {
                type_tag: "DataQuery".to_string(),
                id: String::new(),
                name: String::new(),
                query: query.to_string(),
                connection_id: String::new(),
                connection: DataConnection {
                    id: String::new(),
                    data_source: String::new(),
                    connection_string: None,
                    credentials: None,
                    parameters: HashMap::new(),
                    tags: HashMap::new(),
                },
                params: func_params,
                options: HashMap::new(),
                aux: HashMap::new(),
            },
            options: HashMap::new(),
            parameter_values: param_values,
            aux: HashMap::new(),
            log: None,
        }
    }

    /// A query with one search-pattern parameter, as the platform sends it:
    /// the param carries `options["pattern"] = <type>`, and the already-parsed matcher
    /// lives in `FuncCall.options["patterns"][<param name>]`.
    pub fn make_pattern_call(
        query: &str,
        name: &str,
        pattern_type: &str,
        matcher: serde_json::Value,
    ) -> FuncCall {
        let mut call = make_call(query, vec![(name, "string", serde_json::json!(""))]);
        call.func.params[0].options.insert("pattern".to_string(), pattern_type.to_string());
        call.options.insert("patterns".to_string(), serde_json::json!({ name: matcher }));
        call
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::*;
    use super::*;

    #[test]
    fn test_no_params() {
        let call = make_call("SELECT * FROM t", vec![]);
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t");
        assert!(prepared.values.is_empty());
    }

    #[test]
    fn test_single_int_param() {
        let call = make_call(
            "SELECT * FROM t WHERE id = @id",
            vec![("id", "int", serde_json::json!(42))],
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE id = 42");
    }

    #[test]
    fn test_string_param() {
        let call = make_call(
            "SELECT * FROM t WHERE name = @name",
            vec![("name", "string", serde_json::json!("O'Brien"))],
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE name = 'O''Brien'");
    }

    #[test]
    fn test_multiple_params() {
        let call = make_call(
            "SELECT * FROM t WHERE id = @id AND name = @name",
            vec![
                ("id", "int", serde_json::json!(1)),
                ("name", "string", serde_json::json!("test")),
            ],
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE id = 1 AND name = 'test'");
    }

    #[test]
    fn test_same_param_twice() {
        let call = make_call(
            "SELECT * FROM t WHERE a = @val OR b = @val",
            vec![("val", "int", serde_json::json!(5))],
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE a = 5 OR b = 5");
    }

    #[test]
    fn test_bool_param() {
        let call = make_call(
            "SELECT * FROM t WHERE active = @active",
            vec![("active", "bool", serde_json::json!(true))],
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE active = TRUE");
    }

    #[test]
    fn test_null_param() {
        let call = make_call(
            "SELECT * FROM t WHERE x = @x",
            vec![("x", "string", serde_json::Value::Null)],
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE x = NULL");
    }

    #[test]
    fn test_float_param() {
        let call = make_call(
            "SELECT * FROM t WHERE price > @price",
            vec![("price", "double", serde_json::json!(9.99))],
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE price > 9.99");
    }

    #[test]
    fn test_numeric_pattern_comparison() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @id(id)",
            "id",
            "int",
            serde_json::json!({"op": ">=", "values": [29]}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM mock_data WHERE (id >= 29)");
    }

    #[test]
    fn test_numeric_pattern_in() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @id(id)",
            "id",
            "int",
            serde_json::json!({"op": "in", "values": [29, 30]}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM mock_data WHERE (id in (29, 30))");
    }

    #[test]
    fn test_numeric_pattern_range() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @id(id)",
            "id",
            "int",
            serde_json::json!({"op": "-", "values": [29, 30]}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM mock_data WHERE (id >= 29 AND id <= 30)");
    }

    #[test]
    fn test_numeric_pattern_none_is_always_true() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @id(id)",
            "id",
            "int",
            serde_json::json!({"op": "none", "values": []}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM mock_data WHERE (1 = 1)");
    }

    #[test]
    fn test_string_pattern_contains_lowercases() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @first_name(first_name)",
            "first_name",
            "string",
            serde_json::json!({"op": "contains", "values": ["Z"]}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM mock_data WHERE (LOWER(first_name) LIKE '%z%')");
    }

    #[test]
    fn test_string_pattern_starts_with() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @first_name(first_name)",
            "first_name",
            "string",
            serde_json::json!({"op": "starts with", "values": ["W"]}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM mock_data WHERE (LOWER(first_name) LIKE 'w%')");
    }

    #[test]
    fn test_string_pattern_in_keeps_original_case() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @country(country)",
            "country",
            "string",
            serde_json::json!({"op": "in", "values": ["Poland", "Brazil"]}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM mock_data WHERE (country in ('Poland', 'Brazil'))");
    }

    /// A dialect with no `regex_match` override must fail the query, not drop the filter:
    /// `(1 = 1)` would return the whole unfiltered table as if it were the filtered result.
    /// Matches Java's `getRegexQuery`, which throws `UnsupportedOperationException`.
    #[test]
    fn test_string_pattern_regex_without_dialect_support_errors() {
        let call = make_pattern_call(
            "SELECT * FROM mock_data WHERE @email(email)",
            "email",
            "string",
            serde_json::json!({"op": "regex", "values": ["^\\w+@google.com.au$"]}),
        );
        let error = prepare_query(&call, &DefaultDialect).unwrap_err();
        assert!(error.to_string().contains("REGEXP is not supported"), "got: {}", error);
    }

    #[test]
    fn test_datetime_pattern_range_inclusive() {
        let call = make_pattern_call(
            "SELECT * FROM t WHERE @d(d)",
            "d",
            "datetime",
            serde_json::json!({
                "op": "range",
                "values": ["2021-01-01T00:00:00", "2022-01-01T00:00:00"],
                "include1": true, "include2": true
            }),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(
            prepared.sql,
            "SELECT * FROM t WHERE (d >= '2021-01-01 00:00:00' AND d <= '2022-01-01 00:00:00')"
        );
    }

    #[test]
    fn test_pattern_without_matcher_applies_no_filter() {
        // No `patterns` entry: must not fall back to inlining the raw expression as `'>1'(id)`.
        let mut call = make_call("SELECT * FROM t WHERE @id(id)", vec![("id", "string", serde_json::json!(">1"))]);
        call.func.params[0].options.insert("pattern".to_string(), "int".to_string());
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE (1 = 1)");
    }

    #[test]
    fn test_pattern_and_plain_param_coexist() {
        let mut call = make_call(
            "SELECT * FROM t WHERE @id(id) AND bool = @bool",
            vec![("id", "string", serde_json::json!(">1")), ("bool", "bool", serde_json::json!(false))],
        );
        call.func.params[0].options.insert("pattern".to_string(), "int".to_string());
        call.options.insert(
            "patterns".to_string(),
            serde_json::json!({"id": {"op": ">", "values": [1]}}),
        );
        let prepared = prepare_query(&call, &DefaultDialect).unwrap();
        assert_eq!(prepared.sql, "SELECT * FROM t WHERE (id > 1) AND bool = FALSE");
    }
}
