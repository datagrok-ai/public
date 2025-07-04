import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace queries {
  export async function clickHouseUuidType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseUuidType', {});
  }

  export async function clickHouseSignedIntTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseSignedIntTypes', {});
  }

  export async function clickHouseUnsignedIntTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseUnsignedIntTypes', {});
  }

  export async function clickHouseFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseFloatTypes', {});
  }

  export async function clickHouseArrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseArrayType', {});
  }

  export async function clickHouseTupleType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseTupleType', {});
  }

  export async function clickHouseMapType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseMapType', {});
  }

  export async function clickHouseDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseDateTypes', {});
  }

  export async function clickHouseNestedType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseNestedType', {});
  }

  export async function clickHouseGeoType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseGeoType', {});
  }

  export async function clickHousePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHousePatternsAll', {});
  }

  export async function clickHouseIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseIntTypePatternNone', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpMore', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpMoreEq', { id });
  }

  export async function clickHouseStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpLessEq', { id });
  }

  export async function clickHouseStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpLess', { id });
  }

  export async function clickHouseStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpIn', { id });
  }

  export async function clickHouseStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpNotIn', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpMinMax', { id });
  }

  export async function clickHouseStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypeIntPatternOpNotEq', { id });
  }

  export async function clickHouseDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseDoubleTypePatternNone', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternDoubleOpMore', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternDoubleOpLess', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function clickHouseStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternStringOpContains', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternStringOpIn', { country });
  }

  export async function clickHouseStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHouseStringTypePatternStringOpRegex', { email });
  }

  export async function clickHousePatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:ClickHousePatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function mariaDbDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbDateTypes', {});
  }

  export async function mariaDbBinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbBinaryTypes', {});
  }

  export async function mariaDbBitType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbBitType', {});
  }

  export async function mariaDbCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbCharacterTypes', {});
  }

  export async function mariaDbFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbFloatTypes', {});
  }

  export async function mariaDbIntegerTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbIntegerTypes', {});
  }

  export async function mariaDbSpatialTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbSpatialTypes', {});
  }

  export async function mariaDbPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbPatternsAll', {});
  }

  export async function mariaDbIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbIntTypePatternNone', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpMore', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mariaDbStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpLessEq', { id });
  }

  export async function mariaDbStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpLess', { id });
  }

  export async function mariaDbStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpIn', { id });
  }

  export async function mariaDbStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpNotIn', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpMinMax', { id });
  }

  export async function mariaDbStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypeIntPatternOpNotEq', { id });
  }

  export async function mariaDbDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbDoubleTypePatternNone', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mariaDbStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternStringOpContains', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternStringOpIn', { country });
  }

  export async function mariaDbStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbStringTypePatternStringOpRegex', { email });
  }

  export async function mariaDbPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MariaDbPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function mssqldateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLDateTypes', {});
  }

  export async function mssqlcharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLCharacterTypes', {});
  }

  export async function mssqlxmlType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLXmlType', {});
  }

  export async function mssqlbinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLBinaryTypes', {});
  }

  export async function mssqlnumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLNumericTypes', {});
  }

  export async function mssqlfloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLFloatTypes', {});
  }

  export async function mssqlmoneyTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLMoneyTypes', {});
  }

  export async function mssqlgeographyType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLGeographyType', {});
  }

  export async function mssqlgeometryType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLGeometryType', {});
  }

  export async function mssqlpatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLPatternsAll', {});
  }

  export async function mssqlintTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLIntTypePatternNone', { id });
  }

  export async function mssqlstringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpMore', { id });
  }

  export async function mssqlstringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mssqlstringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpLessEq', { id });
  }

  export async function mssqlstringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpLess', { id });
  }

  export async function mssqlstringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpIn', { id });
  }

  export async function mssqlstringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpNotIn', { id });
  }

  export async function mssqlstringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpMinMax', { id });
  }

  export async function mssqlstringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypeIntPatternOpNotEq', { id });
  }

  export async function mssqldoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLDoubleTypePatternNone', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mssqlstringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternStringOpContains', { first_name });
  }

  export async function mssqlstringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mssqlstringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mssqlstringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLStringTypePatternStringOpIn', { country });
  }

  export async function mssqlpatternsAllParams(first_name: string, id: string, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLPatternsAllParams', { first_name, id, email, some_number, country, date });
  }

  export async function mssqlall(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLAll', {});
  }

  export async function mssqlbyInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLByInt', { orderid });
  }

  export async function mssqlbyStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLByStringPatternInt', { shipVia });
  }

  export async function mssqlbyDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLByDouble', { freight });
  }

  export async function mssqlbyStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLByStringChoices', { shipCountry });
  }

  export async function mssqlbyStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLByStringPatternString', { shipCity });
  }

  export async function mssqlbyDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLByDatetime', { requiredDate });
  }

  export async function mssqlbyStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLByStringPatternDatetime', { orderDate });
  }

  export async function mssqlorders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function mssqlproducts(ProductID: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MSSQLProducts', { ProductID });
  }

  export async function mySqlDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlDateTypes', {});
  }

  export async function mySqlBinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlBinaryTypes', {});
  }

  export async function mySqlBitType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlBitType', {});
  }

  export async function mySqlCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlCharacterTypes', {});
  }

  export async function mySqlFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlFloatTypes', {});
  }

  export async function mySqlJsonType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlJsonType', {});
  }

  export async function mySqlIntegerTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlIntegerTypes', {});
  }

  export async function mySqlSpatialTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlSpatialTypes', {});
  }

  export async function mySqlPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlPatternsAll', {});
  }

  export async function mySqlIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlIntTypePatternNone', { id });
  }

  export async function mySqlStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpMore', { id });
  }

  export async function mySqlStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mySqlStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpLessEq', { id });
  }

  export async function mySqlStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpLess', { id });
  }

  export async function mySqlStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpIn', { id });
  }

  export async function mySqlStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpNotIn', { id });
  }

  export async function mySqlStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpMinMax', { id });
  }

  export async function mySqlStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypeIntPatternOpNotEq', { id });
  }

  export async function mySqlDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlDoubleTypePatternNone', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mySqlStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternStringOpContains', { first_name });
  }

  export async function mySqlStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mySqlStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mySqlStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternStringOpIn', { country });
  }

  export async function mySqlStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlStringTypePatternStringOpRegex', { email });
  }

  export async function mySqlPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:MySqlPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function oracleXmlType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleXmlType', {});
  }

  export async function oracleCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleCharacterTypes', {});
  }

  export async function oracleDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleDateTypes', {});
  }

  export async function oracleJsonType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleJsonType', {});
  }

  export async function oracleUriType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleUriType', {});
  }

  export async function oracleVarrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleVarrayType', {});
  }

  export async function oracleNumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleNumericTypes', {});
  }

  export async function oracleLobTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleLobTypes', {});
  }

  export async function oraclePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OraclePatternsAll', {});
  }

  export async function oracleIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleIntTypePatternNone', { id });
  }

  export async function oracleStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpMore', { id });
  }

  export async function oracleStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpMoreEq', { id });
  }

  export async function oracleStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpLessEq', { id });
  }

  export async function oracleStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpLess', { id });
  }

  export async function oracleStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpIn', { id });
  }

  export async function oracleStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpNotIn', { id });
  }

  export async function oracleStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpMinMax', { id });
  }

  export async function oracleStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypeIntPatternOpNotEq', { id });
  }

  export async function oracleDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleDoubleTypePatternNone', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternDoubleOpMore', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternDoubleOpLess', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function oracleStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternStringOpContains', { first_name });
  }

  export async function oracleStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function oracleStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function oracleStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternStringOpIn', { country });
  }

  export async function oracleStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OracleStringTypePatternStringOpRegex', { email });
  }

  export async function oraclePatternsAllParams(first_name: string, id: string, email: string, some_number: string, country: string, dat: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:OraclePatternsAllParams', { first_name, id, email, some_number, country, dat });
  }

  export async function packageTest(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PackageTest', {});
  }

  export async function postgresAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresAll', {});
  }

  export async function postgresByInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByInt', { orderid });
  }

  export async function postgresByStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByStringPatternInt', { shipVia });
  }

  export async function postgresByDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByDouble', { freight });
  }

  export async function postgresByStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByStringChoices', { shipCountry });
  }

  export async function postgresByStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByStringPatternString', { shipCity });
  }

  export async function postgresByBool(freightLess100: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByBool', { freightLess100 });
  }

  export async function postgresByDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByDatetime', { requiredDate });
  }

  export async function postgresByStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresByStringPatternDatetime', { orderDate });
  }

  export async function postgresOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function postgresProducts(ProductID: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresProducts', { ProductID });
  }

  export async function postgresqlTestCacheTableNormal(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlTestCacheTableNormal', {});
  }

  export async function postgresqlScalarCacheTest(): Promise<number> {
    return await grok.data.query('DBTests:PostgresqlScalarCacheTest', {});
  }

  export async function postgresqlCachedConnTest(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlCachedConnTest', {});
  }

  export async function postgresqlCacheInvalidateQueryTest(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlCacheInvalidateQueryTest', {});
  }

  export async function postgresqlScalarCacheInvalidationTest(): Promise<number> {
    return await grok.data.query('DBTests:PostgresqlScalarCacheInvalidationTest', {});
  }

  export async function postgresqlGetNow(): Promise<any> {
    return await grok.data.query('DBTests:PostgresqlGetNow', {});
  }

  export async function postgresqlScalarCacheTestClient(): Promise<number> {
    return await grok.data.query('DBTests:PostgresqlScalarCacheTestClient', {});
  }

  export async function postgresqlTestCacheTableNormalClient(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlTestCacheTableNormalClient', {});
  }

  export async function testConnCache(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:TestConnCache', {});
  }

  export async function postgresqlTableWideCachedClient(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlTableWideCachedClient', {});
  }

  export async function postgresqlTableWideCachedServer(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlTableWideCachedServer', {});
  }

  export async function postgresqlArrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlArrayType', {});
  }

  export async function postgresqlBitString1(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlBitString1', {});
  }

  export async function postgresqlBitString2(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlBitString2', {});
  }

  export async function postgresqlComposite(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlComposite', {});
  }

  export async function postgresqlDates(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlDates', {});
  }

  export async function postgresqlJSONB(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlJSONB', {});
  }

  export async function postgresqlNumeric(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlNumeric', {});
  }

  export async function postgresqlDouble(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlDouble', {});
  }

  export async function postgresqlReal(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlReal', {});
  }

  export async function postgresqlBigInt(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlBigInt', {});
  }

  export async function postgresqlNumericPrecision(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlNumericPrecision', {});
  }

  export async function postgresqlNumericPrecisionScale(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlNumericPrecisionScale', {});
  }

  export async function postgresqlSmallSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlSmallSerial', {});
  }

  export async function postgresqlSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlSerial', {});
  }

  export async function postgresqlBigSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlBigSerial', {});
  }

  export async function postgresqlUUID(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlUUID', {});
  }

  export async function postgresqlXML(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlXML', {});
  }

  export async function postgresqlPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlPatternsAll', {});
  }

  export async function postgresqlIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlIntTypePatternNone', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpMore', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpMoreEq', { id });
  }

  export async function postgresqlStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpLessEq', { id });
  }

  export async function postgresqlStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpLess', { id });
  }

  export async function postgresqlStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpIn', { id });
  }

  export async function postgresqlStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpNotIn', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpMinMax', { id });
  }

  export async function postgresqlStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypeIntPatternOpNotEq', { id });
  }

  export async function postgresqlDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlDoubleTypePatternNone', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternDoubleOpMore', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternDoubleOpLess', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function postgresqlStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternStringOpContains', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternStringOpIn', { country });
  }

  export async function postgresqlStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlStringTypePatternStringOpRegex', { email });
  }

  export async function postgresqlPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function postgresqlAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlAll', {});
  }

  export async function postgresqlByInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByInt', { orderid });
  }

  export async function postgresqlByStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByStringPatternInt', { shipVia });
  }

  export async function postgresqlByDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByDouble', { freight });
  }

  export async function postgresqlByStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByStringChoices', { shipCountry });
  }

  export async function postgresqlByStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByStringPatternString', { shipCity });
  }

  export async function postgresqlByBool(freightLess100: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByBool', { freightLess100 });
  }

  export async function postgresqlByDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByDatetime', { requiredDate });
  }

  export async function postgresqlByStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlByStringPatternDatetime', { orderDate });
  }

  export async function postgresqlOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function postgresqlProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlProducts', {});
  }

  export async function postgresPerf(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresPerf', {});
  }

  export async function postgreScalarOutput(): Promise<number> {
    return await grok.data.query('DBTests:PostgreScalarOutput', {});
  }

  export async function postgresqlTableNormal(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlTableNormal', {});
  }

  export async function postgresqlTableWide(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlTableWide', {});
  }

  export async function postgresqlTableLong(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:PostgresqlTableLong', {});
  }

  export async function testForColumnsOnEmptyResult(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:TestForColumnsOnEmptyResult', {});
  }

  export async function simpleSelect(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SimpleSelect', {});
  }

  export async function snowflakeBinaryType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeBinaryType', {});
  }

  export async function snowflakeDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeDateTypes', {});
  }

  export async function snowflakeGeoType(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeGeoType', {});
  }

  export async function snowflakeNumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeNumericTypes', {});
  }

  export async function snowflakeSemiStructuredTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeSemiStructuredTypes', {});
  }

  export async function snowflakePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakePatternsAll', {});
  }

  export async function snowflakeIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeIntTypePatternNone', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpMore', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpMoreEq', { id });
  }

  export async function snowflakeStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpLessEq', { id });
  }

  export async function snowflakeStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpLess', { id });
  }

  export async function snowflakeStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpIn', { id });
  }

  export async function snowflakeStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpNotIn', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpMinMax', { id });
  }

  export async function snowflakeStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypeIntPatternOpNotEq', { id });
  }

  export async function snowflakeDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeDoubleTypePatternNone', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternDoubleOpMore', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternDoubleOpLess', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function snowflakeStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternStringOpContains', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternStringOpIn', { country });
  }

  export async function snowflakeStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakeStringTypePatternStringOpRegex', { email });
  }

  export async function snowflakePatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('DBTests:SnowflakePatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }
}

export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('DBTests:Info', {});
  }

  export async function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): Promise<any> {
    return await grok.functions.call('DBTests:ExpectTable', { actual, expected });
  }
}
