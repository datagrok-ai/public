--name: OracleXmlType
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleXmlType(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_xml_type.d42'))
SELECT * FROM xml_data
--end

--name: OracleCharacterTypes
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleCharacterTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_character_types.d42'))
SELECT TRIM(CH) AS CH, VARCH, TRIM(NCH) AS NCH, NVARCH FROM character_type
--end

--name: OracleDateTypes
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleDateTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_date_types.d42'))
SELECT * FROM dates_type
--end

--name: OracleJsonType
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleJsonType(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_json_type.d42'))
SELECT * FROM json_data
--end

--name: OracleUriType
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleUriType(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_uri_type.d42'))
SELECT u.uri.getURL() AS URI FROM uri_types u
--end

--name: OracleVarrayType
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleVarrayType(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_varray_type.d42'))
SELECT * FROM varrays
--end

--name: OracleNumericTypes
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleNumericTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_numeric_types.d42'))
SELECT * FROM numeric_type
--end

--name: OracleLobTypes
--connection: OracleApiTests
--test: ApiTests:expectTable(OracleLobTypes(), OpenFile('System:AppData/ApiTests/datasets/tests/oracle/oracle_lob_types.d42'))
SELECT UTL_RAW.CAST_TO_VARCHAR2(BLOB_TYPE) AS BLOB_TYPE, TO_CHAR(CLOB_TYPE) AS CLOB_TYPE,
       TO_CHAR(NCLOB_TYPE) AS NCLOB_TYPE FROM lobs
--end
