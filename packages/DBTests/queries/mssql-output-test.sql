--name: MSSQLDateTypes
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLDateTypes(), OpenFile('System:AppData/DBTests/mssql/mssql_output_dates.d42'))
SELECT * FROM date_types;
--end

--name: MSSQLCharacterTypes
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLCharacterTypes(), OpenFile('System:AppData/DBTests/mssql/mssql_output_characters.d42'))
SELECT * FROM character_types;
--end

--name: MSSQLXmlType
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLXmlType(), OpenFile('System:AppData/DBTests/mssql/mssql_output_xml.d42'))
SELECT * FROM xml_type;
--end

--name: MSSQLBinaryTypes
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLBinaryTypes(), OpenFile('System:AppData/DBTests/mssql/mssql_output_binary.d42'))
SELECT CONVERT(VARCHAR(max), binary_data, 0) as binary_data, CONVERT(int, varbinary_data, 1)
    as varbinary_data FROM binary_types;
--end

--name: MSSQLNumericTypes
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLNumericTypes(), OpenFile('System:AppData/DBTests/mssql/mssql_output_numeric.d42'))
SELECT * FROM numeric_types;
--end

--name: MSSQLFloatTypes
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLFloatTypes(), OpenFile('System:AppData/DBTests/mssql/mssql_output_float.d42'))
SELECT * FROM float_types;
--end

--name: MSSQLMoneyTypes
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLMoneyTypes(), OpenFile('System:AppData/DBTests/mssql/mssql_output_money.d42'))
SELECT * FROM money_types;
--end

--name: MSSQLGeographyType
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLGeographyType(), OpenFile('System:AppData/DBTests/mssql/mssql_output_geog.d42'))
SELECT id, GeogCol1.ToString() as GeogCol1, GeogCol2 FROM SpatialTable2;
--end

--name: MSSQLGeometryType
--connection: MSSQLDBTests
--test: DBTests:expectTable(MSSQLGeometryType(), OpenFile('System:AppData/DBTests/mssql/mssql_output_geom.d42'))
SELECT id, GeomCol1.ToString() as GeomCol1, GeomCol2 FROM SpatialTable1;
--end
