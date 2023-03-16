--name: RedshiftBinaryTypes
--connection: RedshiftApiTests
--test: ApiTests:expectTable(RedshiftBinaryTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/redshift_output_binary.d42'))
SELECT from_varbyte(varbyte_type, 'utf-8') varbyte_type, from_varbyte(varbinary_type, 'utf-8')
    varbinary_type FROM BINARY_TYPES;
--end

--name: RedshiftCharacterTypes
--connection: RedshiftApiTests
--test: ApiTests:expectTable(RedshiftCharacterTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/redshift_output_characters.d42'))
SELECT * FROM CHARACTER_TYPES;
--end

--name: RedshiftDateTypes
--connection: RedshiftApiTests
--test: ApiTests:expectTable(RedshiftDateTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/redshift_output_dates.d42'))
SELECT * FROM DATE_TYPES;
--end

--name: RedshiftIntegerTypes
--connection: RedshiftApiTests
--test: ApiTests:expectTable(RedshiftIntegerTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/redshift_output_integer.d42'))
SELECT * FROM INTEGER_TYPES;
--end

--name: RedshiftFloatTypes
--connection: RedshiftApiTests
--test: ApiTests:expectTable(RedshiftFloatTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/redshift_output_float.d42'))
SELECT * FROM FLOAT_TYPES;
--end

--name: RedshiftSuperType
--connection: RedshiftApiTests
--test: ApiTests:expectTable(RedshiftSuperType(), OpenFile('System:AppData/ApiTests/datasets/tests/redshift/redshift_output_super.d42'))
SELECT * FROM SUPER_TYPE;
--end
