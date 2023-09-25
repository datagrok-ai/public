--name: MSSQLDateTypes
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLDateTypes(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_dates.d42')) // cat: MSSQL
SELECT * FROM date_types;
--end
``
--name: MSSQLCharacterTypes
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLCharacterTypes(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_characters.d42')) // cat: MSSQL
SELECT * FROM character_types;
--end

--name: MSSQLXmlType
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLXmlType(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_xml.d42')) // cat: MSSQL
SELECT * FROM xml_type;
--end

--name: MSSQLBinaryTypes
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLBinaryTypes(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_binary.d42')) // cat: MSSQL
SELECT CONVERT(VARCHAR(max), binary_data, 0) as binary_data, CONVERT(int, varbinary_data, 1)
    as varbinary_data FROM binary_types;
--end

--name: MSSQLNumericTypes
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLNumericTypes(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_numeric.d42')) // cat: MSSQL
SELECT * FROM numeric_types;
--end

--name: MSSQLFloatTypes
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLFloatTypes(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_float.d42')) // cat: MSSQL
SELECT * FROM float_types;
--end

--name: MSSQLMoneyTypes
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLMoneyTypes(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_money.d42')) // cat: MSSQL
SELECT * FROM money_types;
--end

--name: MSSQLGeographyType
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLGeographyType(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_geog.d42')) // cat: MSSQL
SELECT id, GeogCol1.ToString() as GeogCol1, GeogCol2 FROM SpatialTable2;
--end

--name: MSSQLGeometryType
--connection: MSSQLDBTests
--test: Dbtests:expectTable(MSSQLGeometryType(), OpenFile('System:AppData/Dbtests/mssql/mssql_output_geom.d42')) // cat: MSSQL
SELECT id, GeomCol1.ToString() as GeomCol1, GeomCol2 FROM SpatialTable1;
--end
