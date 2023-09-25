-- name: ClickHouseUuidType
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseUuidType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_uuid.d42')) // cat: ClickHouse
SELECT * FROM uuid_type;
-- end

-- name: ClickHouseSignedIntTypes
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseSignedIntTypes(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_signed_ints.d42')) // cat: ClickHouse
SELECT * FROM SIGNED_INTEGER_TYPES ORDER BY int8_type DESC;
-- end

-- name: ClickHouseUnsignedIntTypes
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseUnsignedIntTypes(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_unsigned_ints.d42')) // cat: ClickHouse
SELECT * FROM UNSIGNED_INTEGER_TYPES;
-- end

-- name: ClickHouseFloatTypes
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseFloatTypes(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_floats.d42')) // cat: ClickHouse
SELECT * FROM FLOAT_TYPES ORDER BY float32_type DESC;
-- end

-- name: ClickHouseArrayType
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseArrayType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_array.d42')) // cat: ClickHouse
SELECT * FROM ARRAY_TYPE;
-- end

-- name: ClickHouseTupleType
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseTupleType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_tuple.d42')) // cat: ClickHouse
SELECT * FROM TUPLE_TYPE;
-- end

-- name: ClickHouseMapType
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseMapType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_map.d42')) // cat: ClickHouse
SELECT * FROM MAP_TYPE;
-- end

-- name: ClickHouseDateTypes
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseDateTypes(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_dates.d42')) // cat: ClickHouse
SELECT * FROM DATE_TYPES;
-- end

-- name: ClickHouseNestedType
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseNestedType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_nested.d42')) // cat: ClickHouse
SELECT * FROM NESTED_TYPE;
-- end

-- name: ClickHouseGeoType
-- connection: ClickHouseDBTests
-- test: Dbtests:expectTable(ClickHouseGeoType(), OpenFile('System:AppData/Dbtests/clickhouse/clickhouse_output_geo.d42')) // cat: ClickHouse
SELECT * FROM GEO_TYPES;
-- end
