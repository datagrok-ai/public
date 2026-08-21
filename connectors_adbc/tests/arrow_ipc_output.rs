use arrow::record_batch::RecordBatch;
use arrow_array::*;
use arrow_schema::{DataType, Field, Schema};
use arrow_ipc::reader::StreamReader;
use std::sync::Arc;

use grok_connect_adbc::arrow_ipc::{record_batch_to_ipc, record_batches_to_ipc};

fn roundtrip(batch: &RecordBatch) -> Vec<RecordBatch> {
    let ipc_bytes = record_batch_to_ipc(batch).unwrap();
    assert!(!ipc_bytes.is_empty(), "IPC bytes should not be empty");
    let cursor = std::io::Cursor::new(ipc_bytes);
    let reader = StreamReader::try_new(cursor, None).unwrap();
    reader.map(|r| r.unwrap()).collect()
}

#[test]
fn test_int_types() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("i8", DataType::Int8, false),
        Field::new("i16", DataType::Int16, false),
        Field::new("i32", DataType::Int32, false),
        Field::new("i64", DataType::Int64, false),
        Field::new("u8", DataType::UInt8, false),
        Field::new("u16", DataType::UInt16, false),
        Field::new("u32", DataType::UInt32, false),
        Field::new("u64", DataType::UInt64, false),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(Int8Array::from(vec![1i8, -1])),
        Arc::new(Int16Array::from(vec![2i16, -2])),
        Arc::new(Int32Array::from(vec![3i32, -3])),
        Arc::new(Int64Array::from(vec![4i64, -4])),
        Arc::new(UInt8Array::from(vec![5u8, 6])),
        Arc::new(UInt16Array::from(vec![7u16, 8])),
        Arc::new(UInt32Array::from(vec![9u32, 10])),
        Arc::new(UInt64Array::from(vec![11u64, 12])),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded[0].num_rows(), 2);
    assert_eq!(decoded[0].num_columns(), 8);
}

#[test]
fn test_float_types() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("f32", DataType::Float32, false),
        Field::new("f64", DataType::Float64, false),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(Float32Array::from(vec![1.5f32, -2.5])),
        Arc::new(Float64Array::from(vec![3.14, -2.71])),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded[0].num_rows(), 2);
}

#[test]
fn test_string_and_bool() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("str", DataType::Utf8, false),
        Field::new("bool", DataType::Boolean, false),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(StringArray::from(vec!["hello", "world", ""])),
        Arc::new(BooleanArray::from(vec![true, false, true])),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded[0].num_rows(), 3);
}

#[test]
fn test_date_and_timestamp() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("date", DataType::Date32, false),
        Field::new("ts", DataType::Timestamp(arrow_schema::TimeUnit::Microsecond, None), false),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(Date32Array::from(vec![19000, 19001])),
        Arc::new(TimestampMicrosecondArray::from(vec![1000000i64, 2000000])),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded[0].num_rows(), 2);
}

#[test]
fn test_nullable_every_type() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("n_i32", DataType::Int32, true),
        Field::new("n_f64", DataType::Float64, true),
        Field::new("n_str", DataType::Utf8, true),
        Field::new("n_bool", DataType::Boolean, true),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(Int32Array::from(vec![Some(1), None, Some(3)])),
        Arc::new(Float64Array::from(vec![None, Some(2.0), None])),
        Arc::new(StringArray::from(vec![Some("a"), None, Some("c")])),
        Arc::new(BooleanArray::from(vec![Some(true), None, Some(false)])),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded[0].num_rows(), 3);
    assert_eq!(decoded[0].column(0).null_count(), 1);
    assert_eq!(decoded[0].column(1).null_count(), 2);
    assert_eq!(decoded[0].column(2).null_count(), 1);
    assert_eq!(decoded[0].column(3).null_count(), 1);
}

#[test]
fn test_zero_rows() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Int32, false),
        Field::new("name", DataType::Utf8, false),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(Int32Array::from(Vec::<i32>::new())),
        Arc::new(StringArray::from(Vec::<&str>::new())),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded.len(), 1);
    assert_eq!(decoded[0].num_rows(), 0);
}

#[test]
fn test_large_utf8() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("large_str", DataType::LargeUtf8, false),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(LargeStringArray::from(vec!["short", &"x".repeat(100_000)])),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded[0].num_rows(), 2);
}

#[test]
fn test_decimal128() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("decimal", DataType::Decimal128(18, 4), false),
    ]));
    let batch = RecordBatch::try_new(schema, vec![
        Arc::new(Decimal128Array::from(vec![12345678i128, -99990000])
            .with_precision_and_scale(18, 4).unwrap()),
    ]).unwrap();

    let decoded = roundtrip(&batch);
    assert_eq!(decoded[0].num_rows(), 2);
}

#[test]
fn test_multi_batch_stream() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Int32, false),
    ]));
    let batch1 = RecordBatch::try_new(
        schema.clone(),
        vec![Arc::new(Int32Array::from(vec![1, 2]))],
    ).unwrap();
    let batch2 = RecordBatch::try_new(
        schema,
        vec![Arc::new(Int32Array::from(vec![3, 4, 5]))],
    ).unwrap();

    let ipc_bytes = record_batches_to_ipc(&[batch1, batch2]).unwrap();
    let cursor = std::io::Cursor::new(ipc_bytes);
    let reader = StreamReader::try_new(cursor, None).unwrap();
    let batches: Vec<_> = reader.map(|r| r.unwrap()).collect();
    assert_eq!(batches.len(), 2);
    assert_eq!(batches[0].num_rows(), 2);
    assert_eq!(batches[1].num_rows(), 3);
}

#[test]
fn test_ipc_schema_preservation() {
    let schema = Arc::new(Schema::new(vec![
        Field::new("n_nationkey", DataType::Int32, false),
        Field::new("n_name", DataType::Utf8, false),
        Field::new("n_regionkey", DataType::Int32, false),
        Field::new("n_comment", DataType::Utf8, true),
    ]));
    let batch = RecordBatch::try_new(schema.clone(), vec![
        Arc::new(Int32Array::from(vec![0])),
        Arc::new(StringArray::from(vec!["ALGERIA"])),
        Arc::new(Int32Array::from(vec![0])),
        Arc::new(StringArray::from(vec![Some("comment")])),
    ]).unwrap();

    let ipc_bytes = record_batch_to_ipc(&batch).unwrap();
    let cursor = std::io::Cursor::new(ipc_bytes);
    let reader = StreamReader::try_new(cursor, None).unwrap();
    let decoded_schema = reader.schema();
    assert_eq!(decoded_schema.fields().len(), 4);
    assert_eq!(decoded_schema.field(0).name(), "n_nationkey");
    assert_eq!(decoded_schema.field(1).name(), "n_name");
    assert_eq!(decoded_schema.field(2).name(), "n_regionkey");
    assert_eq!(decoded_schema.field(3).name(), "n_comment");
}
