-- name: ClickHouseUuidType
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseUuidType(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_uuid.d42'))
SELECT * FROM uuid_type;
-- end

-- name: ClickHouseSignedIntTypes
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseSignedIntTypes(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_signed_ints.d42'))
SELECT * FROM SIGNED_INTEGER_TYPES ORDER BY int8_type DESC;
-- end

-- name: ClickHouseUnsignedIntTypes
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseUnsignedIntTypes(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_unsigned_ints.d42'))
SELECT * FROM UNSIGNED_INTEGER_TYPES;
-- end

-- name: ClickHouseFloatTypes
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseFloatTypes(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_floats.d42'))
SELECT * FROM FLOAT_TYPES ORDER BY float32_type DESC;
-- end

-- name: ClickHouseArrayType
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseArrayType(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_array.d42'))
SELECT * FROM ARRAY_TYPE;
-- end

-- name: ClickHouseTupleType
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseTupleType(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_tuple.d42'))
SELECT * FROM TUPLE_TYPE;
-- end

-- name: ClickHouseMapType
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseMapType(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_map.d42'))
SELECT * FROM MAP_TYPE;
-- end

-- name: ClickHouseDateTypes
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseDateTypes(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_dates.d42'))
SELECT * FROM DATE_TYPES;
-- end

-- name: ClickHouseNestedType
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseNestedType(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_nested.d42'))
SELECT * FROM NESTED_TYPE;
-- end

-- name: ClickHouseGeoType
-- connection: ClickHouseDBTests
-- test: DBTests:expectTable(ClickHouseGeoType(), OpenFile('System:AppData/DBTests/clickhouse/clickhouse_output_geo.d42'))
SELECT * FROM GEO_TYPES;
-- end
