-- Multi-chunk streaming coverage for the ClickHouse ADBC/Arrow connector.
--
-- The Rust connector cuts a chunk every MAX_CHUNK_ROWS = 50_000 rows
-- (grok_connect_adbc/src/routes/query_socket.rs), so every query here returns >= 200_000 rows
-- and exercises the full SIZE -> ACK -> binary -> PART_OK pipeline rather than the
-- single-chunk path the type-coverage suites in clickhouse-adbc-*.sql take.
--
-- All row sources are the `numbers()` generator, not stored tables: the queries are
-- self-contained, need no DDL on the shared demo DB, and stay valid if the `test` database
-- is ever reloaded. Flip the connection's `connectionType` to JDBC to run the identical SQL
-- through the Java connector for an apples-to-apples comparison.
--
-- No `-- test:` annotations on purpose. These have no d42 baselines (a 1M-row baseline is not
-- worth committing, and `generateUUIDv4` is non-deterministic); they are driven from the
-- `ClickHouse Streaming` benchmark category in src/benchmarks/clickhouse-streaming.ts.

-- name: ClickHouseStreamingInts
-- friendlyName: Streaming: 1M ints
-- connection: ClickHouseStreaming
-- Narrowest possible row: isolates chunking/round-trip overhead from payload size. ~20 chunks.
SELECT number AS id FROM numbers(1000000);
-- end

-- name: ClickHouseStreamingLowCardinality
-- friendlyName: Streaming: 1M low-cardinality strings
-- connection: ClickHouseStreaming
-- LowCardinality(String) arrives as an Arrow Dictionary only because the connector sets
-- output_format_arrow_low_cardinality_as_dictionary; this is the query that proves the setting
-- survives across every chunk of a pooled session, not just the first.
SELECT toLowCardinality(['Fail', 'Pending', 'Success'][(number % 3) + 1]) AS status
FROM numbers(1000000);
-- end

-- name: ClickHouseStreamingMixed
-- friendlyName: Streaming: 1M mixed types
-- connection: ClickHouseStreaming
-- The realistic shape: dictionary + plain string + narrow Date/DateTime + float, 1M rows.
-- Date arrives as UInt16 and DateTime as UInt32, so this drives repair_temporal_columns on
-- every chunk — the per-chunk transform most likely to break under streaming.
SELECT
  toUInt32(number) AS id,
  toLowCardinality(['Fail', 'Pending', 'Success'][(number % 3) + 1]) AS status,
  concat('user_', toString(number % 10000)) AS name,
  toDate('2020-01-01') + toIntervalDay(number % 2000) AS created_date,
  toDateTime('2020-01-01 00:00:00') + toIntervalSecond(number) AS created_at,
  toFloat64(number) / 7 AS score
FROM numbers(1000000);
-- end

-- name: ClickHouseStreamingWideStrings
-- friendlyName: Streaming: 250K wide strings
-- connection: ClickHouseStreaming
-- 512 distinct chars per row (~128 MB), so dictionary encoding cannot help and the wire is
-- dominated by plain Utf8. This is the case where externalDataFrameCompress on the ADBC
-- endpoint costs the most, and where Arrow IPC is compared against Java's d42 head-on.
SELECT
  number AS id,
  repeat(hex(MD5(toString(number))), 16) AS payload
FROM numbers(250000);
-- end

-- name: ClickHouseStreamingUnsupportedTypes
-- friendlyName: Streaming: 200K UUID + wide ints
-- connection: ClickHouseStreaming
-- UUID and Int128/UInt256 cannot be serialized to Arrow by ClickHouse, so the connector
-- rewrites them to String via SELECT * REPLACE and re-tags the wide ints as bigint. Both
-- transforms are driven off a single DESCRIBE and must stay aligned across all chunks.
SELECT
  number AS id,
  generateUUIDv4() AS uuid_col,
  toInt128(number) * toInt128(1000000000000000000) AS int128_col,
  toUInt256(number) * toUInt256(1000000000000000000) AS uint256_col
FROM numbers(200000);
-- end

-- name: ClickHouseStreamingRows
-- friendlyName: Streaming: N rows
-- connection: ClickHouseStreaming
-- input: int rows = 1000000
-- Dial the row count by hand when hunting for a chunk-boundary bug. toUInt64 because the
-- connector inlines parameters as SQL literals and the platform renders JSON numbers as
-- doubles (1000000.0).
SELECT
  number AS id,
  toLowCardinality(['Fail', 'Pending', 'Success'][(number % 3) + 1]) AS status
FROM numbers(toUInt64(@rows));
-- end
