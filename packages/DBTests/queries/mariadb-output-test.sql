-- name: MariaDbDateTypes
-- connection: MariaDbDBTests
-- test: Dbtests:expectTable(MariaDbDateTypes(), OpenFile('System:AppData/Dbtests/mysql/mysql_output_dates.d42')) // cat: MariaDb
SELECT * FROM DATE_TYPES;
-- end

-- name: MariaDbBinaryTypes
-- connection: MariaDbDBTests
-- test: Dbtests:expectTable(MariaDbBinaryTypes(), OpenFile('System:AppData/Dbtests/mysql/mysql_output_binary.d42')) // cat: MariaDb
SELECT CAST(binary_type as CHAR(5)) as binary_type, CAST(varbinary_type AS CHAR(8))
                                    as varbinary_type FROM BINARY_TYPES;
-- end

-- name: MariaDbBitType
-- connection: MariaDbDBTests
-- test: Dbtests:expectTable(MariaDbBitType(), OpenFile('System:AppData/Dbtests/mysql/mysql_output_bit.d42')) // cat: MariaDb
SELECT BIN(bit_type1) AS bit_type1, bit_type2 FROM BIT_TYPE;
-- end

-- name: MariaDbCharacterTypes
-- connection: MariaDbDBTests
-- test: Dbtests:expectTable(MariaDbCharacterTypes(), OpenFile('System:AppData/Dbtests/mysql/mysql_output_characters.d42')) // cat: MariaDb
SELECT * FROM CHARACTER_TYPES;
-- end

-- name: MariaDbFloatTypes
-- connection: MariaDbDBTests
-- test: Dbtests:expectTable(MariaDbFloatTypes(), OpenFile('System:AppData/Dbtests/mysql/mysql_output_floats.d42')) // cat: MariaDb
SELECT * FROM FLOAT_TYPES;
-- end

-- name: MariaDbIntegerTypes
-- connection: MariaDbDBTests
-- test: Dbtests:expectTable(MariaDbIntegerTypes(), OpenFile('System:AppData/Dbtests/mysql/mysql_output_integers.d42')) // cat: MariaDb
SELECT * FROM INTEGER_TYPES;
-- end

-- name: MariaDbSpatialTypes
-- connection: MariaDbDBTests
-- test: Dbtests:expectTable(MariaDbSpatialTypes(), OpenFile('System:AppData/Dbtests/mysql/mysql_output_spatial.d42')) // cat: MariaDb
SELECT ST_AsText(geometry_type) AS geometry_type, ST_AsText(point_type) AS point_type FROM GEOMETRY_TYPE;
-- end
