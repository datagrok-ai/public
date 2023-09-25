--name: AthenaArrayIntType
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaArrayIntType(), OpenFile('System:AppData/Dbtests/athena/athena_output_array_int.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM array_type_int;
--end

--name: AthenaArrayStringType
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaArrayStringType(), OpenFile('System:AppData/Dbtests/athena/athena_output_array_string.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM array_type_string;
--end

--name: AthenaCharacterTypes
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaCharacterTypes(), OpenFile('System:AppData/Dbtests/athena/athena_output_characters.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM character_types;
--end

--name: AthenaDateTypes
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaDateTypes(), OpenFile('System:AppData/Dbtests/athena/athena_output_dates.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM date_types;
--end

--name: AthenaFloatTypes
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaFloatTypes(), OpenFile('System:AppData/Dbtests/athena/athena_output_floats.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM float_types;
--end

--name: AthenaMapType
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaMapType(), OpenFile('System:AppData/Dbtests/athena/athena_output_map.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM map_type;
--end

--name: AthenaNumericTypes
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaNumericTypes(), OpenFile('System:AppData/Dbtests/athena/athena_output_numeric.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM numeric_types;
--end

--name: AthenaStructType
--connection: AthenaDBTests
--test: Dbtests:expectTable(AthenaStructType(), OpenFile('System:AppData/Dbtests/athena/athena_output_struct.d42')) //skip:GROK-12423, cat: Athena
SELECT * FROM struct_type;
--end
