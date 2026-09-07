//! Utf8 → `Dictionary<Int32, Utf8>` encoding for low-cardinality string columns.
//!
//! Some ADBC drivers emit plain `Utf8` for every string column regardless of
//! cardinality — they have no equivalent of a `LowCardinality` type to carry the hint,
//! so a repeated value costs its full length on every row. For workloads where the same
//! 10–100k distinct values repeat across 1M rows (`bench_categorical`), that is 2–3×
//! more wire bytes than necessary. Java
//! `grok_connect`'s d42 format dict-encodes strings via `StringColumn.categorize()`;
//! we do the equivalent here before Arrow IPC serialisation so all three downstream
//! legs (Datlas gzip, browser gunzip, `Arrow.fromFeather`) operate on the smaller,
//! already-deduplicated payload. The Datagrok Arrow JS decoder has a fast path for
//! `Dictionary<Int32, Utf8>` at `public/packages/Arrow/src/api/api.ts:199`
//! (`stringColumnFromDictionary`) — no JS changes needed.

use std::collections::HashSet;
use std::sync::Arc;

use anyhow::Result;
use arrow::array::{Array, StringArray};
use arrow::compute::cast;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

#[derive(Debug, Clone, Copy)]
pub struct DictEncodeConfig {
    /// Encode when `n_distinct / n_rows <= max_cardinality_ratio`. Default 0.5.
    pub max_cardinality_ratio: f64,
    /// Skip batches with fewer rows than this. Default 256.
    pub min_rows_for_probe: usize,
    /// Abort the probe as soon as `n_distinct / n_rows > early_exit_ratio`. Default 0.5.
    pub early_exit_ratio: f64,
}

impl Default for DictEncodeConfig {
    fn default() -> Self {
        Self {
            max_cardinality_ratio: 0.5,
            min_rows_for_probe: 256,
            early_exit_ratio: 0.5,
        }
    }
}

impl DictEncodeConfig {
    /// Parse a per-query `meta.dictEncodeStrings` / env `DICT_ENCODE_STRINGS` value:
    ///   - `"off"` / `"false"` / `"0"` → `None` (disable)
    ///   - `"auto"` / `"on"` / `"true"` → default config
    ///   - `"<float>"` (e.g. `"0.3"`) → default with `max_cardinality_ratio = <float>`
    pub fn from_str_option(value: &str) -> Option<Self> {
        let v = value.trim().to_ascii_lowercase();
        match v.as_str() {
            "off" | "false" | "0" | "no" => None,
            "" | "auto" | "on" | "true" | "yes" => Some(Self::default()),
            other => match other.parse::<f64>() {
                Ok(ratio) if ratio > 0.0 && ratio <= 1.0 => Some(Self {
                    max_cardinality_ratio: ratio,
                    early_exit_ratio: ratio,
                    ..Self::default()
                }),
                _ => Some(Self::default()),
            },
        }
    }
}

/// Returns a new `RecordBatch` with qualifying `Utf8` columns rebuilt as
/// `Dictionary<Int32, Utf8>`. Non-qualifying columns pass through unchanged
/// (same `ArrayRef`, no copy).
///
/// Qualification per column:
/// 1. Must be `Utf8`.
/// 2. `num_rows >= config.min_rows_for_probe`.
/// 3. Cardinality probe (short-circuited by `early_exit_ratio`) finds
///    `n_distinct / n_rows <= config.max_cardinality_ratio`.
/// 4. At least one non-null value (all-null columns are left untouched —
///    there is no benefit to encoding them and it would only grow the schema).
pub fn maybe_encode_batch(batch: &RecordBatch, config: DictEncodeConfig) -> Result<RecordBatch> {
    let n_rows = batch.num_rows();
    let n_cols = batch.num_columns();

    let mut new_columns = Vec::with_capacity(n_cols);
    let mut new_fields = Vec::with_capacity(n_cols);
    let mut any_encoded = false;
    let mut bytes_before: usize = 0;
    let mut bytes_after: usize = 0;

    let schema = batch.schema();
    for (idx, field) in schema.fields().iter().enumerate() {
        let col = batch.column(idx);
        let col_bytes_before = col.get_array_memory_size();
        bytes_before += col_bytes_before;

        let should_encode = matches!(field.data_type(), DataType::Utf8)
            && n_rows >= config.min_rows_for_probe
            && {
                let s = col.as_any().downcast_ref::<StringArray>().expect("Utf8");
                probe_worth_encoding(s, config)
            };

        if should_encode {
            let dict_type = DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8));
            let encoded = cast(col, &dict_type)?;
            let col_bytes_after = encoded.get_array_memory_size();
            bytes_after += col_bytes_after;
            let new_field = Arc::new(Field::new(
                field.name(),
                dict_type,
                field.is_nullable(),
            ).with_metadata(field.metadata().clone()));
            tracing::debug!(
                column = %field.name(),
                n_rows,
                bytes_before = col_bytes_before,
                bytes_after = col_bytes_after,
                "Dict-encoded Utf8 column"
            );
            new_fields.push(new_field);
            new_columns.push(encoded);
            any_encoded = true;
        }
        else {
            bytes_after += col_bytes_before;
            new_fields.push(field.clone());
            new_columns.push(col.clone());
        }
    }

    if !any_encoded {
        return Ok(batch.clone());
    }

    tracing::info!(
        n_rows,
        bytes_before,
        bytes_after,
        "Dict encoding applied to chunk"
    );

    let new_schema = Arc::new(Schema::new_with_metadata(
        new_fields,
        schema.metadata().clone(),
    ));
    Ok(RecordBatch::try_new(new_schema, new_columns)?)
}

/// Fast cardinality probe with early exit. Returns `true` iff the column has
/// at least one non-null value *and* the distinct-value ratio among non-null
/// rows is `<= config.max_cardinality_ratio`.
fn probe_worth_encoding(s: &StringArray, config: DictEncodeConfig) -> bool {
    let n_rows = s.len();
    if n_rows == 0 {
        return false;
    }
    // Distinct-count cap we're willing to tolerate. Compare against the full
    // row count so the threshold is row-count-based and stable regardless of
    // how many nulls appear — matches the intuitive "n_distinct / n_rows".
    let max_allowed = ((n_rows as f64) * config.max_cardinality_ratio).floor() as usize;
    let early_exit_cap = ((n_rows as f64) * config.early_exit_ratio).floor() as usize;
    let cap = max_allowed.max(early_exit_cap);

    let mut seen: HashSet<&str> = HashSet::with_capacity(cap.min(1024).max(16));
    for i in 0..n_rows {
        if s.is_null(i) {
            continue;
        }
        let v = s.value(i);
        if seen.insert(v) && seen.len() > cap {
            return false;
        }
    }
    if seen.is_empty() {
        return false;
    }
    seen.len() <= max_allowed
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{ArrayRef, Float64Array, Int32Array};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow_ipc::reader::StreamReader;
    use std::sync::Arc;
    use std::time::Instant;

    fn utf8_schema(nullable: bool) -> Arc<Schema> {
        Arc::new(Schema::new(vec![Field::new("s", DataType::Utf8, nullable)]))
    }

    fn make_utf8_batch(values: Vec<Option<&str>>) -> RecordBatch {
        let schema = utf8_schema(true);
        let array: ArrayRef = Arc::new(StringArray::from(values));
        RecordBatch::try_new(schema, vec![array]).unwrap()
    }

    #[test]
    fn test_low_cardinality_utf8_gets_encoded() {
        let values: Vec<Option<&str>> = (0..1000)
            .map(|i| Some(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"][i % 10]))
            .collect();
        let batch = make_utf8_batch(values);
        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        assert!(matches!(
            out.schema().field(0).data_type(),
            DataType::Dictionary(_, _)
        ));
        // Round-trip back to Utf8 to compare values.
        let back = cast(out.column(0), &DataType::Utf8).unwrap();
        let back = back.as_any().downcast_ref::<StringArray>().unwrap();
        for i in 0..1000 {
            assert_eq!(back.value(i), ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"][i % 10]);
        }
    }

    #[test]
    fn test_high_cardinality_utf8_unchanged() {
        let values: Vec<String> = (0..1000).map(|i| format!("uuid-{:08x}", i)).collect();
        let refs: Vec<Option<&str>> = values.iter().map(|s| Some(s.as_str())).collect();
        let batch = make_utf8_batch(refs);
        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        assert_eq!(*out.schema().field(0).data_type(), DataType::Utf8);
    }

    #[test]
    fn test_mixed_schema_only_affects_utf8() {
        let low_card: Vec<Option<&str>> = (0..1000)
            .map(|i| Some(["red", "green", "blue"][i % 3]))
            .collect();
        let schema = Arc::new(Schema::new(vec![
            Field::new("i", DataType::Int32, false),
            Field::new("f", DataType::Float64, false),
            Field::new("s", DataType::Utf8, false),
        ]));
        let int_arr: ArrayRef = Arc::new(Int32Array::from((0..1000).collect::<Vec<_>>()));
        let float_arr: ArrayRef = Arc::new(Float64Array::from((0..1000).map(|i| i as f64).collect::<Vec<_>>()));
        let str_arr: ArrayRef = Arc::new(StringArray::from(low_card));
        let batch = RecordBatch::try_new(schema, vec![int_arr.clone(), float_arr.clone(), str_arr]).unwrap();

        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        // Int32 and Float64 columns should pass through as the *same* ArrayRef.
        assert!(Arc::ptr_eq(out.column(0), &int_arr), "int column should be Arc::ptr_eq");
        assert!(Arc::ptr_eq(out.column(1), &float_arr), "float column should be Arc::ptr_eq");
        // Utf8 column should be rebuilt as a dictionary.
        assert!(matches!(out.schema().field(2).data_type(), DataType::Dictionary(_, _)));
    }

    #[test]
    fn test_nullable_strings_preserve_null_count() {
        let values: Vec<Option<&str>> = (0..1000)
            .map(|i| if i % 3 == 0 { None } else { Some(["x", "y", "z", "w", "v"][i % 5]) })
            .collect();
        let nulls_expected = values.iter().filter(|v| v.is_none()).count();
        let batch = make_utf8_batch(values.clone());
        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        assert!(matches!(out.schema().field(0).data_type(), DataType::Dictionary(_, _)));
        assert_eq!(out.column(0).null_count(), nulls_expected);

        let back = cast(out.column(0), &DataType::Utf8).unwrap();
        let back = back.as_any().downcast_ref::<StringArray>().unwrap();
        for (i, expected) in values.iter().enumerate() {
            match expected {
                None => assert!(back.is_null(i), "row {} should be null", i),
                Some(s) => assert_eq!(back.value(i), *s),
            }
        }
    }

    #[test]
    fn test_empty_batch() {
        let schema = utf8_schema(true);
        let array: ArrayRef = Arc::new(StringArray::from(Vec::<Option<&str>>::new()));
        let batch = RecordBatch::try_new(schema.clone(), vec![array]).unwrap();
        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        assert_eq!(out.num_rows(), 0);
        assert_eq!(*out.schema().field(0).data_type(), DataType::Utf8);
    }

    #[test]
    fn test_small_batch_below_min_rows() {
        let values: Vec<Option<&str>> = (0..100).map(|_| Some("only-one")).collect();
        let batch = make_utf8_batch(values);
        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        // 100 < min_rows_for_probe (256) → pass-through.
        assert_eq!(*out.schema().field(0).data_type(), DataType::Utf8);
    }

    #[test]
    fn test_threshold_boundary() {
        // 1000 rows, exactly 500 distinct → ratio 0.50 → encoded (<=).
        let vals_eq: Vec<String> = (0..1000).map(|i| format!("v{}", i / 2)).collect();
        let refs_eq: Vec<Option<&str>> = vals_eq.iter().map(|s| Some(s.as_str())).collect();
        let batch_eq = make_utf8_batch(refs_eq);
        let out_eq = maybe_encode_batch(&batch_eq, DictEncodeConfig::default()).unwrap();
        assert!(
            matches!(out_eq.schema().field(0).data_type(), DataType::Dictionary(_, _)),
            "500/1000 = 0.50 should encode (ratio boundary is inclusive)"
        );

        // 1000 rows, 510 distinct → ratio 0.51 → not encoded.
        let vals_gt: Vec<String> = (0..1000)
            .map(|i| {
                if i < 980 { format!("v{}", i / 2) } else { format!("uniq{}", i) }
            })
            .collect();
        let refs_gt: Vec<Option<&str>> = vals_gt.iter().map(|s| Some(s.as_str())).collect();
        let batch_gt = make_utf8_batch(refs_gt);
        let out_gt = maybe_encode_batch(&batch_gt, DictEncodeConfig::default()).unwrap();
        assert_eq!(*out_gt.schema().field(0).data_type(), DataType::Utf8);
    }

    #[test]
    fn test_all_null_utf8() {
        let values: Vec<Option<&str>> = (0..1000).map(|_| None).collect();
        let batch = make_utf8_batch(values);
        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        assert_eq!(*out.schema().field(0).data_type(), DataType::Utf8);
    }

    #[test]
    fn test_round_trip_via_ipc() {
        let values: Vec<Option<&str>> = (0..1000)
            .map(|i| Some(["alpha", "beta", "gamma"][i % 3]))
            .collect();
        let batch = make_utf8_batch(values.clone());
        let encoded = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        assert!(matches!(
            encoded.schema().field(0).data_type(),
            DataType::Dictionary(_, _)
        ));

        let ipc = crate::arrow_ipc::record_batches_to_ipc(&[encoded]).unwrap();
        let cursor = std::io::Cursor::new(ipc);
        let mut reader = StreamReader::try_new(cursor, None).unwrap();
        let decoded = reader.next().unwrap().unwrap();
        assert!(matches!(
            decoded.schema().field(0).data_type(),
            DataType::Dictionary(_, _)
        ));
        let back = cast(decoded.column(0), &DataType::Utf8).unwrap();
        let back = back.as_any().downcast_ref::<StringArray>().unwrap();
        for (i, expected) in values.iter().enumerate() {
            assert_eq!(back.value(i), expected.unwrap());
        }
    }

    #[test]
    fn test_early_exit_keeps_probe_cheap() {
        let values: Vec<String> = (0..1_000_000).map(|i| format!("uniq-{:010}", i)).collect();
        let refs: Vec<Option<&str>> = values.iter().map(|s| Some(s.as_str())).collect();
        let batch = make_utf8_batch(refs);
        let start = Instant::now();
        let out = maybe_encode_batch(&batch, DictEncodeConfig::default()).unwrap();
        let elapsed = start.elapsed();
        assert_eq!(*out.schema().field(0).data_type(), DataType::Utf8);
        // Early exit should fire after ~500k distinct, well before the full 1M scan.
        // Give generous headroom for CI variance; the naïve path would be >> 500ms.
        assert!(
            elapsed.as_millis() < 500,
            "probe took {:?}, expected < 500ms thanks to early exit",
            elapsed
        );
    }
}
