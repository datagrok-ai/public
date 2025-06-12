import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Queries {
  export async function clickHouseUuidType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseUuidType', {});
  }

  export async function clickHouseSignedIntTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseSignedIntTypes', {});
  }

  export async function clickHouseUnsignedIntTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseUnsignedIntTypes', {});
  }

  export async function clickHouseFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseFloatTypes', {});
  }

  export async function clickHouseArrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseArrayType', {});
  }

  export async function clickHouseTupleType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseTupleType', {});
  }

  export async function clickHouseMapType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseMapType', {});
  }

  export async function clickHouseDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseDateTypes', {});
  }

  export async function clickHouseNestedType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseNestedType', {});
  }

  export async function clickHouseGeoType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseGeoType', {});
  }

  export async function clickHousePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHousePatternsAll', {});
  }

  export async function clickHouseIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseIntTypePatternNone', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpMore', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpMoreEq', { id });
  }

  export async function clickHouseStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpLessEq', { id });
  }

  export async function clickHouseStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpLess', { id });
  }

  export async function clickHouseStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpIn', { id });
  }

  export async function clickHouseStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpNotIn', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpMinMax', { id });
  }

  export async function clickHouseStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypeIntPatternOpNotEq', { id });
  }

  export async function clickHouseDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseDoubleTypePatternNone', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternDoubleOpMore', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternDoubleOpLess', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function clickHouseStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternStringOpContains', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternStringOpIn', { country });
  }

  export async function clickHouseStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHouseStringTypePatternStringOpRegex', { email });
  }

  export async function clickHousePatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:ClickHousePatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function mariaDbDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbDateTypes', {});
  }

  export async function mariaDbBinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbBinaryTypes', {});
  }

  export async function mariaDbBitType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbBitType', {});
  }

  export async function mariaDbCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbCharacterTypes', {});
  }

  export async function mariaDbFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbFloatTypes', {});
  }

  export async function mariaDbIntegerTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbIntegerTypes', {});
  }

  export async function mariaDbSpatialTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbSpatialTypes', {});
  }

  export async function mariaDbPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbPatternsAll', {});
  }

  export async function mariaDbIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbIntTypePatternNone', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpMore', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mariaDbStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpLessEq', { id });
  }

  export async function mariaDbStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpLess', { id });
  }

  export async function mariaDbStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpIn', { id });
  }

  export async function mariaDbStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpNotIn', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpMinMax', { id });
  }

  export async function mariaDbStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypeIntPatternOpNotEq', { id });
  }

  export async function mariaDbDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbDoubleTypePatternNone', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mariaDbStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternStringOpContains', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternStringOpIn', { country });
  }

  export async function mariaDbStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbStringTypePatternStringOpRegex', { email });
  }

  export async function mariaDbPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MariaDbPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function mssqldateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLDateTypes', {});
  }

  export async function mssqlcharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLCharacterTypes', {});
  }

  export async function mssqlxmlType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLXmlType', {});
  }

  export async function mssqlbinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLBinaryTypes', {});
  }

  export async function mssqlnumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLNumericTypes', {});
  }

  export async function mssqlfloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLFloatTypes', {});
  }

  export async function mssqlmoneyTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLMoneyTypes', {});
  }

  export async function mssqlgeographyType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLGeographyType', {});
  }

  export async function mssqlgeometryType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLGeometryType', {});
  }

  export async function mssqlpatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLPatternsAll', {});
  }

  export async function mssqlintTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLIntTypePatternNone', { id });
  }

  export async function mssqlstringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpMore', { id });
  }

  export async function mssqlstringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mssqlstringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpLessEq', { id });
  }

  export async function mssqlstringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpLess', { id });
  }

  export async function mssqlstringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpIn', { id });
  }

  export async function mssqlstringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpNotIn', { id });
  }

  export async function mssqlstringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpMinMax', { id });
  }

  export async function mssqlstringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypeIntPatternOpNotEq', { id });
  }

  export async function mssqldoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLDoubleTypePatternNone', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mssqlstringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternStringOpContains', { first_name });
  }

  export async function mssqlstringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mssqlstringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mssqlstringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLStringTypePatternStringOpIn', { country });
  }

  export async function mssqlpatternsAllParams(first_name: string, id: string, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLPatternsAllParams', { first_name, id, email, some_number, country, date });
  }

  export async function mssqlall(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLAll', {});
  }

  export async function mssqlbyInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLByInt', { orderid });
  }

  export async function mssqlbyStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLByStringPatternInt', { shipVia });
  }

  export async function mssqlbyDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLByDouble', { freight });
  }

  export async function mssqlbyStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLByStringChoices', { shipCountry });
  }

  export async function mssqlbyStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLByStringPatternString', { shipCity });
  }

  export async function mssqlbyDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLByDatetime', { requiredDate });
  }

  export async function mssqlbyStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLByStringPatternDatetime', { orderDate });
  }

  export async function mssqlorders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function mssqlproducts(ProductID: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MSSQLProducts', { ProductID });
  }

  export async function mySqlDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlDateTypes', {});
  }

  export async function mySqlBinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlBinaryTypes', {});
  }

  export async function mySqlBitType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlBitType', {});
  }

  export async function mySqlCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlCharacterTypes', {});
  }

  export async function mySqlFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlFloatTypes', {});
  }

  export async function mySqlJsonType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlJsonType', {});
  }

  export async function mySqlIntegerTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlIntegerTypes', {});
  }

  export async function mySqlSpatialTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlSpatialTypes', {});
  }

  export async function mySqlPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlPatternsAll', {});
  }

  export async function mySqlIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlIntTypePatternNone', { id });
  }

  export async function mySqlStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpMore', { id });
  }

  export async function mySqlStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mySqlStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpLessEq', { id });
  }

  export async function mySqlStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpLess', { id });
  }

  export async function mySqlStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpIn', { id });
  }

  export async function mySqlStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpNotIn', { id });
  }

  export async function mySqlStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpMinMax', { id });
  }

  export async function mySqlStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypeIntPatternOpNotEq', { id });
  }

  export async function mySqlDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlDoubleTypePatternNone', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mySqlStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternStringOpContains', { first_name });
  }

  export async function mySqlStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mySqlStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mySqlStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternStringOpIn', { country });
  }

  export async function mySqlStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlStringTypePatternStringOpRegex', { email });
  }

  export async function mySqlPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:MySqlPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function oracleXmlType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleXmlType', {});
  }

  export async function oracleCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleCharacterTypes', {});
  }

  export async function oracleDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleDateTypes', {});
  }

  export async function oracleJsonType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleJsonType', {});
  }

  export async function oracleUriType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleUriType', {});
  }

  export async function oracleVarrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleVarrayType', {});
  }

  export async function oracleNumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleNumericTypes', {});
  }

  export async function oracleLobTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleLobTypes', {});
  }

  export async function oraclePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OraclePatternsAll', {});
  }

  export async function oracleIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleIntTypePatternNone', { id });
  }

  export async function oracleStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpMore', { id });
  }

  export async function oracleStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpMoreEq', { id });
  }

  export async function oracleStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpLessEq', { id });
  }

  export async function oracleStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpLess', { id });
  }

  export async function oracleStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpIn', { id });
  }

  export async function oracleStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpNotIn', { id });
  }

  export async function oracleStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpMinMax', { id });
  }

  export async function oracleStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypeIntPatternOpNotEq', { id });
  }

  export async function oracleDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleDoubleTypePatternNone', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternDoubleOpMore', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternDoubleOpLess', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function oracleStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternStringOpContains', { first_name });
  }

  export async function oracleStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function oracleStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function oracleStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternStringOpIn', { country });
  }

  export async function oracleStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OracleStringTypePatternStringOpRegex', { email });
  }

  export async function oraclePatternsAllParams(first_name: string, id: string, email: string, some_number: string, country: string, dat: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:OraclePatternsAllParams', { first_name, id, email, some_number, country, dat });
  }

  export async function packageTest(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PackageTest', {});
  }

  export async function postgresAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresAll', {});
  }

  export async function postgresByInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByInt', { orderid });
  }

  export async function postgresByStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByStringPatternInt', { shipVia });
  }

  export async function postgresByDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByDouble', { freight });
  }

  export async function postgresByStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByStringChoices', { shipCountry });
  }

  export async function postgresByStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByStringPatternString', { shipCity });
  }

  export async function postgresByBool(freightLess100: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByBool', { freightLess100 });
  }

  export async function postgresByDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByDatetime', { requiredDate });
  }

  export async function postgresByStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresByStringPatternDatetime', { orderDate });
  }

  export async function postgresOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function postgresProducts(ProductID: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresProducts', { ProductID });
  }

  export async function postgresqlTestCacheTableNormal(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlTestCacheTableNormal', {});
  }

  export async function postgresqlScalarCacheTest(): Promise<number> {
    return await grok.data.query('Dbtests:PostgresqlScalarCacheTest', {});
  }

  export async function postgresqlCachedConnTest(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlCachedConnTest', {});
  }

  export async function postgresqlCacheInvalidateQueryTest(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlCacheInvalidateQueryTest', {});
  }

  export async function postgresqlScalarCacheInvalidationTest(): Promise<number> {
    return await grok.data.query('Dbtests:PostgresqlScalarCacheInvalidationTest', {});
  }

  export async function postgresqlGetNow(): Promise<any> {
    return await grok.data.query('Dbtests:PostgresqlGetNow', {});
  }

  export async function postgresqlScalarCacheTestClient(): Promise<number> {
    return await grok.data.query('Dbtests:PostgresqlScalarCacheTestClient', {});
  }

  export async function postgresqlTestCacheTableNormalClient(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlTestCacheTableNormalClient', {});
  }

  export async function testConnCache(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:TestConnCache', {});
  }

  export async function postgresqlTableWideCachedClient(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlTableWideCachedClient', {});
  }

  export async function postgresqlTableWideCachedServer(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlTableWideCachedServer', {});
  }

  export async function postgresqlArrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlArrayType', {});
  }

  export async function postgresqlBitString1(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlBitString1', {});
  }

  export async function postgresqlBitString2(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlBitString2', {});
  }

  export async function postgresqlComposite(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlComposite', {});
  }

  export async function postgresqlDates(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlDates', {});
  }

  export async function postgresqlJSONB(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlJSONB', {});
  }

  export async function postgresqlNumeric(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlNumeric', {});
  }

  export async function postgresqlDouble(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlDouble', {});
  }

  export async function postgresqlReal(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlReal', {});
  }

  export async function postgresqlBigInt(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlBigInt', {});
  }

  export async function postgresqlNumericPrecision(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlNumericPrecision', {});
  }

  export async function postgresqlNumericPrecisionScale(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlNumericPrecisionScale', {});
  }

  export async function postgresqlSmallSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlSmallSerial', {});
  }

  export async function postgresqlSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlSerial', {});
  }

  export async function postgresqlBigSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlBigSerial', {});
  }

  export async function postgresqlUUID(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlUUID', {});
  }

  export async function postgresqlXML(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlXML', {});
  }

  export async function postgresqlPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlPatternsAll', {});
  }

  export async function postgresqlIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlIntTypePatternNone', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpMore', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpMoreEq', { id });
  }

  export async function postgresqlStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpLessEq', { id });
  }

  export async function postgresqlStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpLess', { id });
  }

  export async function postgresqlStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpIn', { id });
  }

  export async function postgresqlStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpNotIn', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpMinMax', { id });
  }

  export async function postgresqlStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypeIntPatternOpNotEq', { id });
  }

  export async function postgresqlDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlDoubleTypePatternNone', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternDoubleOpMore', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternDoubleOpLess', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function postgresqlStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternStringOpContains', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternStringOpIn', { country });
  }

  export async function postgresqlStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlStringTypePatternStringOpRegex', { email });
  }

  export async function postgresqlPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function postgresqlAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlAll', {});
  }

  export async function postgresqlByInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByInt', { orderid });
  }

  export async function postgresqlByStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByStringPatternInt', { shipVia });
  }

  export async function postgresqlByDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByDouble', { freight });
  }

  export async function postgresqlByStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByStringChoices', { shipCountry });
  }

  export async function postgresqlByStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByStringPatternString', { shipCity });
  }

  export async function postgresqlByBool(freightLess100: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByBool', { freightLess100 });
  }

  export async function postgresqlByDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByDatetime', { requiredDate });
  }

  export async function postgresqlByStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlByStringPatternDatetime', { orderDate });
  }

  export async function postgresqlOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function postgresqlProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlProducts', {});
  }

  export async function postgresPerf(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresPerf', {});
  }

  export async function postgreScalarOutput(): Promise<number> {
    return await grok.data.query('Dbtests:PostgreScalarOutput', {});
  }

  export async function postgresqlTableNormal(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlTableNormal', {});
  }

  export async function postgresqlTableWide(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlTableWide', {});
  }

  export async function postgresqlTableLong(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:PostgresqlTableLong', {});
  }

  export async function testForColumnsOnEmptyResult(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:TestForColumnsOnEmptyResult', {});
  }

  export async function simpleSelect(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SimpleSelect', {});
  }

  export async function snowflakeBinaryType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeBinaryType', {});
  }

  export async function snowflakeDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeDateTypes', {});
  }

  export async function snowflakeGeoType(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeGeoType', {});
  }

  export async function snowflakeNumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeNumericTypes', {});
  }

  export async function snowflakeSemiStructuredTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeSemiStructuredTypes', {});
  }

  export async function snowflakePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakePatternsAll', {});
  }

  export async function snowflakeIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeIntTypePatternNone', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpMore', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpMoreEq', { id });
  }

  export async function snowflakeStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpLessEq', { id });
  }

  export async function snowflakeStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpLess', { id });
  }

  export async function snowflakeStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpIn', { id });
  }

  export async function snowflakeStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpNotIn', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpMinMax', { id });
  }

  export async function snowflakeStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypeIntPatternOpNotEq', { id });
  }

  export async function snowflakeDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeDoubleTypePatternNone', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternDoubleOpMore', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternDoubleOpLess', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function snowflakeStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternStringOpContains', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternStringOpIn', { country });
  }

  export async function snowflakeStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakeStringTypePatternStringOpRegex', { email });
  }

  export async function snowflakePatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('Dbtests:SnowflakePatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }
}

export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Dbtests:Info', {});
  }

  export async function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): Promise<any> {
    return await grok.functions.call('Dbtests:ExpectTable', { actual, expected });
  }
}
