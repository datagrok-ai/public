--name: OracleXmlType
--connection: OracleDBTests
--test: DBTests:expectTable(OracleXmlType(), OpenFile('System:AppData/DBTests/oracle/oracle_xml_type.d42'))
SELECT * FROM xml_data
--end

--name: OracleCharacterTypes
--connection: OracleDBTests
--test: DBTests:expectTable(OracleCharacterTypes(), OpenFile('System:AppData/DBTests/oracle/oracle_character_types.d42'))
SELECT TRIM(CH) AS CH, VARCH, TRIM(NCH) AS NCH, NVARCH FROM character_type
--end

--name: OracleDateTypes
--connection: OracleDBTests
--test: DBTests:expectTable(OracleDateTypes(), OpenFile('System:AppData/DBTests/oracle/oracle_date_types.d42'))
SELECT * FROM dates_type
--end

--name: OracleJsonType
--connection: OracleDBTests
--test: DBTests:expectTable(OracleJsonType(), OpenFile('System:AppData/DBTests/oracle/oracle_json_type.d42'))
SELECT * FROM json_data
--end

--name: OracleUriType
--connection: OracleDBTests
--test: DBTests:expectTable(OracleUriType(), OpenFile('System:AppData/DBTests/oracle/oracle_uri_type.d42'))
SELECT u.uri.getURL() AS URI FROM uri_types u
--end

--name: OracleVarrayType
--connection: OracleDBTests
--test: DBTests:expectTable(OracleVarrayType(), OpenFile('System:AppData/DBTests/oracle/oracle_varray_type.d42'))
SELECT * FROM varrays
--end

--name: OracleNumericTypes
--connection: OracleDBTests
--test: DBTests:expectTable(OracleNumericTypes(), OpenFile('System:AppData/DBTests/oracle/oracle_numeric_types.d42'))
SELECT * FROM numeric_type
--end

--name: OracleLobTypes
--connection: OracleDBTests
--test: DBTests:expectTable(OracleLobTypes(), OpenFile('System:AppData/DBTests/oracle/oracle_lob_types.d42'))
SELECT UTL_RAW.CAST_TO_VARCHAR2(BLOB_TYPE) AS BLOB_TYPE, TO_CHAR(CLOB_TYPE) AS CLOB_TYPE,
       TO_CHAR(NCLOB_TYPE) AS NCLOB_TYPE FROM lobs
--end
