--name: AthenaArrayIntType
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaArrayIntType(), OpenFile('System:AppData/DBTests/athena/athena_output_array_int.d42')) //skip:GROK-12423
SELECT * FROM array_type_int;
--end

--name: AthenaArrayStringType
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaArrayStringType(), OpenFile('System:AppData/DBTests/athena/athena_output_array_string.d42')) //skip:GROK-12423
SELECT * FROM array_type_string;
--end

--name: AthenaCharacterTypes
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaCharacterTypes(), OpenFile('System:AppData/DBTests/athena/athena_output_characters.d42')) //skip:GROK-12423
SELECT * FROM character_types;
--end

--name: AthenaDateTypes
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaDateTypes(), OpenFile('System:AppData/DBTests/athena/athena_output_dates.d42')) //skip:GROK-12423
SELECT * FROM date_types;
--end

--name: AthenaFloatTypes
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaFloatTypes(), OpenFile('System:AppData/DBTests/athena/athena_output_floats.d42')) //skip:GROK-12423
SELECT * FROM float_types;
--end

--name: AthenaMapType
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaMapType(), OpenFile('System:AppData/DBTests/athena/athena_output_map.d42')) //skip:GROK-12423
SELECT * FROM map_type;
--end

--name: AthenaNumericTypes
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaNumericTypes(), OpenFile('System:AppData/DBTests/athena/athena_output_numeric.d42')) //skip:GROK-12423
SELECT * FROM numeric_types;
--end

--name: AthenaStructType
--connection: AthenaDBTests
--test: DBTests:expectTable(AthenaStructType(), OpenFile('System:AppData/DBTests/athena/athena_output_struct.d42')) //skip:GROK-12423
SELECT * FROM struct_type;
--end
