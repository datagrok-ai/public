use arrow::record_batch::RecordBatch;
use arrow_ipc::writer::StreamWriter;
use anyhow::Result;

/// Serialize a RecordBatch to Arrow IPC Stream format bytes.
/// The client's Arrow.fromFeather() uses arrow.tableFromIPC() which accepts
/// both Stream and File format. Stream format is preferred — simpler, no footer.
pub fn record_batch_to_ipc(batch: &RecordBatch) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    let mut writer = StreamWriter::try_new(&mut buf, &batch.schema())?;
    writer.write(batch)?;
    writer.finish()?;
    Ok(buf)
}

/// Serialize multiple RecordBatches into a single Arrow IPC stream.
/// Used when accumulating small batches into one chunk.
pub fn record_batches_to_ipc(batches: &[RecordBatch]) -> Result<Vec<u8>> {
    if batches.is_empty() {
        return Ok(Vec::new());
    }
    let schema = batches[0].schema();
    let mut buf = Vec::new();
    let mut writer = StreamWriter::try_new(&mut buf, &schema)?;
    for batch in batches {
        writer.write(batch)?;
    }
    writer.finish()?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_array::{Int32Array, StringArray, Float64Array, BooleanArray, RecordBatch};
    use arrow_schema::{Schema, Field, DataType};
    use arrow_ipc::reader::StreamReader;
    use std::sync::Arc;

    #[test]
    fn test_basic_ipc_roundtrip() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("id", DataType::Int32, false),
            Field::new("name", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(Int32Array::from(vec![1, 2, 3])),
                Arc::new(StringArray::from(vec!["a", "b", "c"])),
            ],
        ).unwrap();

        let ipc_bytes = record_batch_to_ipc(&batch).unwrap();
        assert!(!ipc_bytes.is_empty());

        // Verify it's readable
        let cursor = std::io::Cursor::new(ipc_bytes);
        let mut reader = StreamReader::try_new(cursor, None).unwrap();
        let decoded = reader.next().unwrap().unwrap();
        assert_eq!(decoded.num_rows(), 3);
        assert_eq!(decoded.num_columns(), 2);
    }

    #[test]
    fn test_empty_batch() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("id", DataType::Int32, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![Arc::new(Int32Array::from(Vec::<i32>::new()))],
        ).unwrap();

        let ipc_bytes = record_batch_to_ipc(&batch).unwrap();
        assert!(!ipc_bytes.is_empty());

        let cursor = std::io::Cursor::new(ipc_bytes);
        let mut reader = StreamReader::try_new(cursor, None).unwrap();
        let decoded = reader.next().unwrap().unwrap();
        assert_eq!(decoded.num_rows(), 0);
    }

    #[test]
    fn test_nullable_columns() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("n_int", DataType::Int32, true),
            Field::new("n_str", DataType::Utf8, true),
            Field::new("n_dbl", DataType::Float64, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(Int32Array::from(vec![Some(1), None, Some(3)])),
                Arc::new(StringArray::from(vec![Some("a"), None, Some("c")])),
                Arc::new(Float64Array::from(vec![Some(1.1), Some(2.2), None])),
            ],
        ).unwrap();

        let ipc_bytes = record_batch_to_ipc(&batch).unwrap();
        let cursor = std::io::Cursor::new(ipc_bytes);
        let mut reader = StreamReader::try_new(cursor, None).unwrap();
        let decoded = reader.next().unwrap().unwrap();
        assert_eq!(decoded.num_rows(), 3);
        assert_eq!(decoded.column(0).null_count(), 1);
        assert_eq!(decoded.column(1).null_count(), 1);
        assert_eq!(decoded.column(2).null_count(), 1);
    }

    #[test]
    fn test_all_types() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("i32", DataType::Int32, false),
            Field::new("f64", DataType::Float64, false),
            Field::new("str", DataType::Utf8, false),
            Field::new("bool", DataType::Boolean, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(Int32Array::from(vec![10, 20])),
                Arc::new(Float64Array::from(vec![1.5, 2.5])),
                Arc::new(StringArray::from(vec!["hello", "world"])),
                Arc::new(BooleanArray::from(vec![true, false])),
            ],
        ).unwrap();

        let ipc_bytes = record_batch_to_ipc(&batch).unwrap();
        let cursor = std::io::Cursor::new(ipc_bytes);
        let mut reader = StreamReader::try_new(cursor, None).unwrap();
        let decoded = reader.next().unwrap().unwrap();
        assert_eq!(decoded.num_rows(), 2);
        assert_eq!(decoded.num_columns(), 4);
    }

    #[test]
    fn test_multi_batch_ipc() {
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
        let mut reader = StreamReader::try_new(cursor, None).unwrap();
        let b1 = reader.next().unwrap().unwrap();
        assert_eq!(b1.num_rows(), 2);
        let b2 = reader.next().unwrap().unwrap();
        assert_eq!(b2.num_rows(), 3);
    }
}
