-- Output-type coverage for the experimental ClickHouseADBC (Arrow IPC) connector.
-- Deliberately reuses the SAME baselines captured from the Java ClickHouse connector
-- (System:AppData/Dbtests/clickhouse/clickhouse_output_*.d42): the point is to verify the
-- Arrow -> DataFrame path yields the same DataFrame as the Java d42 -> DataFrame path on
-- identical ClickHouse data. Type mismatches here are the signal, not noise.
--
-- One deliberate exception: FloatTypes targets an ADBC-specific baseline
-- (clickhouse-adbc/clickhouse_output_floats.d42). The Java JDBC path narrows Float64 and every
-- Decimal through a float32, losing precision (e.g. 999.9999 -> 999.9998779296875, and 1.79E+308
-- -> Infinity); the Arrow path preserves full float64. That divergence is ADBC being *more*
-- accurate, not a bug to hide, so it gets its own baseline. Any other type divergence must be
-- fixed in the decode path, not by re-baselining.

-- name: ClickHouseADBCUuidType
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCUuidType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_uuid.d42')) // cat: ClickHouseADBC
SELECT * FROM uuid_type;
-- end

-- name: ClickHouseADBCSignedIntTypes
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCSignedIntTypes(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_signed_ints.d42')) // cat: ClickHouseADBC
SELECT * FROM SIGNED_INTEGER_TYPES ORDER BY int8_type DESC;
-- end

-- name: ClickHouseADBCUnsignedIntTypes
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCUnsignedIntTypes(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_unsigned_ints.d42')) // cat: ClickHouseADBC
SELECT * FROM UNSIGNED_INTEGER_TYPES;
-- end

-- name: ClickHouseADBCFloatTypes
-- connection: ClickHouseADBCDBTests
-- ADBC-specific baseline on purpose: the Arrow path keeps full float64 precision that Java's
-- float32 narrowing loses. See the header note.
-- test: Dbtests:expectTable(ClickHouseADBCFloatTypes(), OpenFile('System:AppData/Dbtests/clickhouse-adbc/clickhouse_output_floats.d42')) // cat: ClickHouseADBC
SELECT * FROM FLOAT_TYPES ORDER BY float32_type DESC;
-- end

-- name: ClickHouseADBCArrayType
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCArrayType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_array.d42')) // cat: ClickHouseADBC
SELECT * FROM ARRAY_TYPE;
-- end

-- name: ClickHouseADBCTupleType
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCTupleType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_tuple.d42')) // cat: ClickHouseADBC
SELECT * FROM TUPLE_TYPE;
-- end

-- name: ClickHouseADBCMapType
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCMapType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_map.d42')) // cat: ClickHouseADBC
SELECT * FROM MAP_TYPE;
-- end

-- name: ClickHouseADBCDateTypes
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCDateTypes(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_dates.d42')) // cat: ClickHouseADBC
SELECT * FROM DATE_TYPES;
-- end

-- name: ClickHouseADBCNestedType
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCNestedType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_nested.d42')) // cat: ClickHouseADBC
SELECT * FROM NESTED_TYPE;
-- end

-- name: ClickHouseADBCGeoType
-- connection: ClickHouseADBCDBTests
-- test: Dbtests:expectTable(ClickHouseADBCGeoType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_geo.d42')) // cat: ClickHouseADBC
SELECT * FROM GEO_TYPES;
-- end
