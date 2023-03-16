--name: MySqlDateTypes
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlDateTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_dates.d42'))
SELECT * FROM DATE_TYPES;
-- end

--name: MySqlBinaryTypes
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlBinaryTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_binary.d42'))
SELECT CAST(binary_type as CHAR(5)) as binary_type, CAST(varbinary_type AS CHAR(8))
    as varbinary_type FROM BINARY_TYPES;
-- end

--name: MySqlBitType
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlBitType(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_bit.d42'))
SELECT BIN(bit_type1) AS bit_type1, bit_type2 FROM BIT_TYPE;
-- end

--name: MySqlCharacterTypes
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlCharacterTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_characters.d42'))
SELECT * FROM CHARACTER_TYPES;
-- end

--name: MySqlFloatTypes
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlFloatTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_floats.d42'))
SELECT * FROM FLOAT_TYPES;
-- end

--name: MySqlJsonType
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlJsonType(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_json.d42'))
SELECT * FROM JSON_TYPE;
-- end

--name: MySqlIntegerTypes
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlIntegerTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_integers.d42'))
SELECT * FROM INTEGER_TYPES;
-- end

--name: MySqlSpatialTypes
--connection: MySqlApiTests
--test: ApiTests:expectTable(MySqlSpatialTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/mysql/mysql_output_spatial.d42'))
SELECT * FROM GEOMETRY_TYPE;
-- end
