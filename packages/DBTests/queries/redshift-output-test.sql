--name: RedshiftBinaryTypes
--connection: RedshiftDBTests
--test: DBTests:expectTable(RedshiftBinaryTypes(), OpenFile('System:AppData/DBTests/redshift/redshift_output_binary.d42')) //skip: GROK-12507
SELECT from_varbyte(varbyte_type, 'utf-8') varbyte_type, from_varbyte(varbinary_type, 'utf-8')
    varbinary_type FROM BINARY_TYPES;
--end

--name: RedshiftCharacterTypes
--connection: RedshiftDBTests
--test: DBTests:expectTable(RedshiftCharacterTypes(), OpenFile('System:AppData/DBTests/redshift/redshift_output_characters.d42')) //skip: GROK-12507
SELECT * FROM CHARACTER_TYPES;
--end

--name: RedshiftDateTypes
--connection: RedshiftDBTests
--test: DBTests:expectTable(RedshiftDateTypes(), OpenFile('System:AppData/DBTests/redshift/redshift_output_dates.d42')) //skip: GROK-12507
SELECT * FROM DATE_TYPES;
--end

--name: RedshiftIntegerTypes
--connection: RedshiftDBTests
--test: DBTests:expectTable(RedshiftIntegerTypes(), OpenFile('System:AppData/DBTests/redshift/redshift_output_integer.d42')) //skip: GROK-12507
SELECT * FROM INTEGER_TYPES;
--end

--name: RedshiftFloatTypes
--connection: RedshiftDBTests
--test: DBTests:expectTable(RedshiftFloatTypes(), OpenFile('System:AppData/DBTests/redshift/redshift_output_float.d42')) //skip: GROK-12507
SELECT * FROM FLOAT_TYPES;
--end

--name: RedshiftSuperType
--connection: RedshiftDBTests
--test: DBTests:expectTable(RedshiftSuperType(), OpenFile('System:AppData/DBTests/redshift/redshift_output_super.d42')) //skip: GROK-12507
SELECT * FROM SUPER_TYPE;
--end
