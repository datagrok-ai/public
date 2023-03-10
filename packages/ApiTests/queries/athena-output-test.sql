--name: AthenaArrayIntType
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaArrayIntType(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_array_int.d42'))
SELECT * FROM array_type_int;
--end

--name: AthenaArrayStringType
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaArrayStringType(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_array_string.d42'))
SELECT * FROM array_type_string;
--end

--name: AthenaCharacterTypes
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaCharacterTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_characters.d42'))
SELECT * FROM character_types;
--end

--name: AthenaDateTypes
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaDateTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_dates.d42'))
SELECT * FROM date_types;
--end

--name: AthenaFloatTypes
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaFloatTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_floats.d42'))
SELECT * FROM float_types;
--end

--name: AthenaMapType
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaMapType(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_map.d42'))
SELECT * FROM map_type;
--end

--name: AthenaNumericTypes
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaNumericTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_numeric.d42'))
SELECT * FROM numeric_types;
--end

--name: AthenaStructType
--connection: AthenaApiTests
--test: ApiTests:expectTable(AthenaStructType(), OpenFile('System:AppData/ApiTests/datasets/tests/athena/athena_output_struct.d42'))
SELECT * FROM struct_type;
--end
