-- name: MariaDbDateTypes
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbDateTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_dates.d42'))
SELECT * FROM DATE_TYPES;
-- end

-- name: MariaDbBinaryTypes
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbBinaryTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_binary.d42'))
SELECT CAST(binary_type as CHAR(5)) as binary_type, CAST(varbinary_type AS CHAR(8))
                                    as varbinary_type FROM BINARY_TYPES;
-- end

-- name: MariaDbBitType
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbBitType(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_bit.d42'))
SELECT BIN(bit_type1) AS bit_type1, bit_type2 FROM BIT_TYPE;
-- end

-- name: MariaDbCharacterTypes
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbCharacterTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_characters.d42'))
SELECT * FROM CHARACTER_TYPES;
-- end

-- name: MariaDbFloatTypes
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbFloatTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_floats.d42'))
SELECT * FROM FLOAT_TYPES;
-- end

-- name: MariaDbIntegerTypes
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbIntegerTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_integers.d42'))
SELECT * FROM INTEGER_TYPES;
-- end

-- name: MariaDbSpatialTypes
-- connection: MariaDbApiTests
-- test: ApiTests:expectTable(MariaDbSpatialTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_spatial.d42'))
SELECT ST_AsText(geometry_type) AS geometry_type, ST_AsText(point_type) AS point_type FROM GEOMETRY_TYPE;
-- end
