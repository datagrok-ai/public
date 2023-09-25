--name: OracleXmlType
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleXmlType(), OpenFile('System:AppData/Dbtests/oracle/oracle_xml_type.d42')) // cat: Oracle
SELECT * FROM xml_data
--end

--name: OracleCharacterTypes
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleCharacterTypes(), OpenFile('System:AppData/Dbtests/oracle/oracle_character_types.d42')) // cat: Oracle
SELECT TRIM(CH) AS CH, VARCH, TRIM(NCH) AS NCH, NVARCH FROM character_type
--end

--name: OracleDateTypes
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleDateTypes(), OpenFile('System:AppData/Dbtests/oracle/oracle_date_types.d42')) // cat: Oracle
SELECT * FROM dates_type
--end

--name: OracleJsonType
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleJsonType(), OpenFile('System:AppData/Dbtests/oracle/oracle_json_type.d42')) // cat: Oracle
SELECT * FROM json_data
--end

--name: OracleUriType
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleUriType(), OpenFile('System:AppData/Dbtests/oracle/oracle_uri_type.d42')) // cat: Oracle
SELECT u.uri.getURL() AS URI FROM uri_types u
--end

--name: OracleVarrayType
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleVarrayType(), OpenFile('System:AppData/Dbtests/oracle/oracle_varray_type.d42')) // cat: Oracle
SELECT * FROM varrays
--end

--name: OracleNumericTypes
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleNumericTypes(), OpenFile('System:AppData/Dbtests/oracle/oracle_numeric_types.d42')) // cat: Oracle
SELECT * FROM numeric_type
--end

--name: OracleLobTypes
--connection: OracleDBTests
--test: Dbtests:expectTable(OracleLobTypes(), OpenFile('System:AppData/Dbtests/oracle/oracle_lob_types.d42')) // cat: Oracle
SELECT UTL_RAW.CAST_TO_VARCHAR2(BLOB_TYPE) AS BLOB_TYPE, TO_CHAR(CLOB_TYPE) AS CLOB_TYPE,
       TO_CHAR(NCLOB_TYPE) AS NCLOB_TYPE FROM lobs
--end
