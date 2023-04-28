-- name: ClickHouseUuidType
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseUuidType(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_uuid.d42'))
SELECT * FROM uuid_type;
-- end

-- name: ClickHouseSignedIntTypes
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseSignedIntTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_signed_ints.d42'))
SELECT * FROM SIGNED_INTEGER_TYPES ORDER BY int8_type DESC;
-- end

-- name: ClickHouseUnsignedIntTypes
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseUnsignedIntTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_unsigned_ints.d42'))
SELECT * FROM UNSIGNED_INTEGER_TYPES;
-- end

-- name: ClickHouseFloatTypes
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseFloatTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_floats.d42'))
SELECT * FROM FLOAT_TYPES ORDER BY float32_type DESC;
-- end

-- name: ClickHouseArrayType
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseArrayType(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_array.d42'))
SELECT * FROM ARRAY_TYPE;
-- end

-- name: ClickHouseTupleType
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseTupleType(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_tuple.d42'))
SELECT * FROM TUPLE_TYPE;
-- end

-- name: ClickHouseMapType
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseMapType(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_map.d42'))
SELECT * FROM MAP_TYPE;
-- end

-- name: ClickHouseDateTypes
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseDateTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_dates.d42'))
SELECT * FROM DATE_TYPES;
-- end

-- name: ClickHouseNestedType
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseNestedType(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_nested.d42'))
SELECT * FROM NESTED_TYPE;
-- end

-- name: ClickHouseGeoType
-- connection: ClickHouseApiTests
-- test: ApiTests:expectTable(ClickHouseGeoType(), OpenFile('System:AppData/ApiTests/datasets/tests/clickhouse/clickhouse_output_geo.d42'))
SELECT * FROM GEO_TYPES;
-- end
