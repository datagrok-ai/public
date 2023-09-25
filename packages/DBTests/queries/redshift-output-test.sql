--name: RedshiftBinaryTypes
--connection: RedshiftDBTests
--test: Dbtests:expectTable(RedshiftBinaryTypes(), OpenFile('System:AppData/Dbtests/redshift/redshift_output_binary.d42')) //skip: GROK-12507, cat: Redshift
SELECT from_varbyte(varbyte_type, 'utf-8') varbyte_type, from_varbyte(varbinary_type, 'utf-8')
    varbinary_type FROM BINARY_TYPES;
--end

--name: RedshiftCharacterTypes
--connection: RedshiftDBTests
--test: Dbtests:expectTable(RedshiftCharacterTypes(), OpenFile('System:AppData/Dbtests/redshift/redshift_output_characters.d42')) //skip: GROK-12507, cat: Redshift
SELECT * FROM CHARACTER_TYPES;
--end

--name: RedshiftDateTypes
--connection: RedshiftDBTests
--test: Dbtests:expectTable(RedshiftDateTypes(), OpenFile('System:AppData/Dbtests/redshift/redshift_output_dates.d42')) //skip: GROK-12507, cat: Redshift
SELECT * FROM DATE_TYPES;
--end

--name: RedshiftIntegerTypes
--connection: RedshiftDBTests
--test: Dbtests:expectTable(RedshiftIntegerTypes(), OpenFile('System:AppData/Dbtests/redshift/redshift_output_integer.d42')) //skip: GROK-12507, cat: Redshift
SELECT * FROM INTEGER_TYPES;
--end

--name: RedshiftFloatTypes
--connection: RedshiftDBTests
--test: Dbtests:expectTable(RedshiftFloatTypes(), OpenFile('System:AppData/Dbtests/redshift/redshift_output_float.d42')) //skip: GROK-12507, cat: Redshift
SELECT * FROM FLOAT_TYPES;
--end

--name: RedshiftSuperType
--connection: RedshiftDBTests
--test: Dbtests:expectTable(RedshiftSuperType(), OpenFile('System:AppData/Dbtests/redshift/redshift_output_super.d42')) //skip: GROK-12507, cat: Redshift
SELECT * FROM SUPER_TYPE;
--end
