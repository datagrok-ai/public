import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace queries {
  export async function clickHouseUuidType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseUuidType', {});
  }

  export async function clickHouseSignedIntTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseSignedIntTypes', {});
  }

  export async function clickHouseUnsignedIntTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseUnsignedIntTypes', {});
  }

  export async function clickHouseFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseFloatTypes', {});
  }

  export async function clickHouseArrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseArrayType', {});
  }

  export async function clickHouseTupleType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseTupleType', {});
  }

  export async function clickHouseMapType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseMapType', {});
  }

  export async function clickHouseDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseDateTypes', {});
  }

  export async function clickHouseNestedType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseNestedType', {});
  }

  export async function clickHouseGeoType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseGeoType', {});
  }

  export async function clickHousePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHousePatternsAll', {});
  }

  export async function clickHouseIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseIntTypePatternNone', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpMore', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpMoreEq', { id });
  }

  export async function clickHouseStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpLessEq', { id });
  }

  export async function clickHouseStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpLess', { id });
  }

  export async function clickHouseStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpIn', { id });
  }

  export async function clickHouseStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpNotIn', { id });
  }

  export async function clickHouseStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpMinMax', { id });
  }

  export async function clickHouseStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypeIntPatternOpNotEq', { id });
  }

  export async function clickHouseDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseDoubleTypePatternNone', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternDoubleOpMore', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternDoubleOpLess', { some_number });
  }

  export async function clickHouseStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function clickHouseStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternStringOpContains', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function clickHouseStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternStringOpIn', { country });
  }

  export async function clickHouseStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHouseStringTypePatternStringOpRegex', { email });
  }

  export async function clickHousePatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:ClickHousePatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function mariaDbDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbDateTypes', {});
  }

  export async function mariaDbBinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbBinaryTypes', {});
  }

  export async function mariaDbBitType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbBitType', {});
  }

  export async function mariaDbCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbCharacterTypes', {});
  }

  export async function mariaDbFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbFloatTypes', {});
  }

  export async function mariaDbIntegerTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbIntegerTypes', {});
  }

  export async function mariaDbSpatialTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbSpatialTypes', {});
  }

  export async function mariaDbPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbPatternsAll', {});
  }

  export async function mariaDbIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbIntTypePatternNone', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpMore', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mariaDbStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpLessEq', { id });
  }

  export async function mariaDbStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpLess', { id });
  }

  export async function mariaDbStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpIn', { id });
  }

  export async function mariaDbStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpNotIn', { id });
  }

  export async function mariaDbStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpMinMax', { id });
  }

  export async function mariaDbStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypeIntPatternOpNotEq', { id });
  }

  export async function mariaDbDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbDoubleTypePatternNone', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mariaDbStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mariaDbStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternStringOpContains', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mariaDbStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternStringOpIn', { country });
  }

  export async function mariaDbStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbStringTypePatternStringOpRegex', { email });
  }

  export async function mariaDbPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MariaDbPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function mssqldateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLDateTypes', {});
  }

  export async function mssqlcharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLCharacterTypes', {});
  }

  export async function mssqlxmlType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLXmlType', {});
  }

  export async function mssqlbinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLBinaryTypes', {});
  }

  export async function mssqlnumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLNumericTypes', {});
  }

  export async function mssqlfloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLFloatTypes', {});
  }

  export async function mssqlmoneyTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLMoneyTypes', {});
  }

  export async function mssqlgeographyType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLGeographyType', {});
  }

  export async function mssqlgeometryType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLGeometryType', {});
  }

  export async function mssqlpatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLPatternsAll', {});
  }

  export async function mssqlintTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLIntTypePatternNone', { id });
  }

  export async function mssqlstringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpMore', { id });
  }

  export async function mssqlstringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mssqlstringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpLessEq', { id });
  }

  export async function mssqlstringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpLess', { id });
  }

  export async function mssqlstringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpIn', { id });
  }

  export async function mssqlstringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpNotIn', { id });
  }

  export async function mssqlstringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpMinMax', { id });
  }

  export async function mssqlstringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypeIntPatternOpNotEq', { id });
  }

  export async function mssqldoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLDoubleTypePatternNone', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mssqlstringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mssqlstringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternStringOpContains', { first_name });
  }

  export async function mssqlstringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mssqlstringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mssqlstringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLStringTypePatternStringOpIn', { country });
  }

  export async function mssqlpatternsAllParams(first_name: string, id: string, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLPatternsAllParams', { first_name, id, email, some_number, country, date });
  }

  export async function mssqlall(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLAll', {});
  }

  export async function mssqlbyInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLByInt', { orderid });
  }

  export async function mssqlbyStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLByStringPatternInt', { shipVia });
  }

  export async function mssqlbyDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLByDouble', { freight });
  }

  export async function mssqlbyStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLByStringChoices', { shipCountry });
  }

  export async function mssqlbyStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLByStringPatternString', { shipCity });
  }

  export async function mssqlbyDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLByDatetime', { requiredDate });
  }

  export async function mssqlbyStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLByStringPatternDatetime', { orderDate });
  }

  export async function mssqlorders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function mssqlproducts(ProductID: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MSSQLProducts', { ProductID });
  }

  export async function mySqlDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlDateTypes', {});
  }

  export async function mySqlBinaryTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlBinaryTypes', {});
  }

  export async function mySqlBitType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlBitType', {});
  }

  export async function mySqlCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlCharacterTypes', {});
  }

  export async function mySqlFloatTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlFloatTypes', {});
  }

  export async function mySqlJsonType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlJsonType', {});
  }

  export async function mySqlIntegerTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlIntegerTypes', {});
  }

  export async function mySqlSpatialTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlSpatialTypes', {});
  }

  export async function mySqlPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlPatternsAll', {});
  }

  export async function mySqlIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlIntTypePatternNone', { id });
  }

  export async function mySqlStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpMore', { id });
  }

  export async function mySqlStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpMoreEq', { id });
  }

  export async function mySqlStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpLessEq', { id });
  }

  export async function mySqlStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpLess', { id });
  }

  export async function mySqlStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpIn', { id });
  }

  export async function mySqlStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpNotIn', { id });
  }

  export async function mySqlStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpMinMax', { id });
  }

  export async function mySqlStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypeIntPatternOpNotEq', { id });
  }

  export async function mySqlDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlDoubleTypePatternNone', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternDoubleOpMore', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternDoubleOpLess', { some_number });
  }

  export async function mySqlStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function mySqlStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternStringOpContains', { first_name });
  }

  export async function mySqlStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function mySqlStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function mySqlStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternStringOpIn', { country });
  }

  export async function mySqlStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlStringTypePatternStringOpRegex', { email });
  }

  export async function mySqlPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:MySqlPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function oracleXmlType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleXmlType', {});
  }

  export async function oracleCharacterTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleCharacterTypes', {});
  }

  export async function oracleDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleDateTypes', {});
  }

  export async function oracleJsonType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleJsonType', {});
  }

  export async function oracleUriType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleUriType', {});
  }

  export async function oracleVarrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleVarrayType', {});
  }

  export async function oracleNumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleNumericTypes', {});
  }

  export async function oracleLobTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleLobTypes', {});
  }

  export async function oraclePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OraclePatternsAll', {});
  }

  export async function oracleIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleIntTypePatternNone', { id });
  }

  export async function oracleStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpMore', { id });
  }

  export async function oracleStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpMoreEq', { id });
  }

  export async function oracleStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpLessEq', { id });
  }

  export async function oracleStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpLess', { id });
  }

  export async function oracleStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpIn', { id });
  }

  export async function oracleStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpNotIn', { id });
  }

  export async function oracleStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpMinMax', { id });
  }

  export async function oracleStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypeIntPatternOpNotEq', { id });
  }

  export async function oracleDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleDoubleTypePatternNone', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternDoubleOpMore', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternDoubleOpLess', { some_number });
  }

  export async function oracleStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function oracleStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternStringOpContains', { first_name });
  }

  export async function oracleStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function oracleStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function oracleStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternStringOpIn', { country });
  }

  export async function oracleStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OracleStringTypePatternStringOpRegex', { email });
  }

  export async function oraclePatternsAllParams(first_name: string, id: string, email: string, some_number: string, country: string, dat: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:OraclePatternsAllParams', { first_name, id, email, some_number, country, dat });
  }

  export async function packageTest(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PackageTest', {});
  }

  export async function postgresAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresAll', {});
  }

  export async function postgresByInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByInt', { orderid });
  }

  export async function postgresByStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByStringPatternInt', { shipVia });
  }

  export async function postgresByDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByDouble', { freight });
  }

  export async function postgresByStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByStringChoices', { shipCountry });
  }

  export async function postgresByStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByStringPatternString', { shipCity });
  }

  export async function postgresByBool(freightLess100: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByBool', { freightLess100 });
  }

  export async function postgresByDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByDatetime', { requiredDate });
  }

  export async function postgresByStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresByStringPatternDatetime', { orderDate });
  }

  export async function postgresOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function postgresProducts(ProductID: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresProducts', { ProductID });
  }

  export async function postgresqlTestCacheTableNormal(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlTestCacheTableNormal', {});
  }

  export async function postgresqlScalarCacheTest(): Promise<number> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlScalarCacheTest', {});
  }

  export async function postgresqlCachedConnTest(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlCachedConnTest', {});
  }

  export async function postgresqlCacheInvalidateQueryTest(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlCacheInvalidateQueryTest', {});
  }

  export async function postgresqlScalarCacheInvalidationTest(): Promise<number> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlScalarCacheInvalidationTest', {});
  }

  export async function postgresqlGetNow(): Promise<any> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlGetNow', {});
  }

  export async function postgresqlScalarCacheTestClient(): Promise<number> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlScalarCacheTestClient', {});
  }

  export async function postgresqlTestCacheTableNormalClient(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlTestCacheTableNormalClient', {});
  }

  export async function testConnCache(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:TestConnCache', {});
  }

  export async function postgresqlTableWideCachedClient(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlTableWideCachedClient', {});
  }

  export async function postgresqlTableWideCachedServer(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlTableWideCachedServer', {});
  }

  export async function postgresqlArrayType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlArrayType', {});
  }

  export async function postgresqlBitString1(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlBitString1', {});
  }

  export async function postgresqlBitString2(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlBitString2', {});
  }

  export async function postgresqlComposite(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlComposite', {});
  }

  export async function postgresqlDates(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlDates', {});
  }

  export async function postgresqlJSONB(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlJSONB', {});
  }

  export async function postgresqlNumeric(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlNumeric', {});
  }

  export async function postgresqlDouble(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlDouble', {});
  }

  export async function postgresqlReal(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlReal', {});
  }

  export async function postgresqlBigInt(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlBigInt', {});
  }

  export async function postgresqlNumericPrecision(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlNumericPrecision', {});
  }

  export async function postgresqlNumericPrecisionScale(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlNumericPrecisionScale', {});
  }

  export async function postgresqlSmallSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlSmallSerial', {});
  }

  export async function postgresqlSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlSerial', {});
  }

  export async function postgresqlBigSerial(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlBigSerial', {});
  }

  export async function postgresqlUUID(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlUUID', {});
  }

  export async function postgresqlXML(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlXML', {});
  }

  export async function postgresqlPatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlPatternsAll', {});
  }

  export async function postgresqlIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlIntTypePatternNone', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpMore', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpMoreEq', { id });
  }

  export async function postgresqlStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpLessEq', { id });
  }

  export async function postgresqlStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpLess', { id });
  }

  export async function postgresqlStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpIn', { id });
  }

  export async function postgresqlStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpNotIn', { id });
  }

  export async function postgresqlStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpMinMax', { id });
  }

  export async function postgresqlStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypeIntPatternOpNotEq', { id });
  }

  export async function postgresqlDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlDoubleTypePatternNone', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternDoubleOpMore', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternDoubleOpLess', { some_number });
  }

  export async function postgresqlStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function postgresqlStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternStringOpContains', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function postgresqlStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternStringOpIn', { country });
  }

  export async function postgresqlStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlStringTypePatternStringOpRegex', { email });
  }

  export async function postgresqlPatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlPatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }

  export async function postgresqlAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlAll', {});
  }

  export async function postgresqlByInt(orderid: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByInt', { orderid });
  }

  export async function postgresqlByStringPatternInt(shipVia: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByStringPatternInt', { shipVia });
  }

  export async function postgresqlByDouble(freight: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByDouble', { freight });
  }

  export async function postgresqlByStringChoices(shipCountry: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByStringChoices', { shipCountry });
  }

  export async function postgresqlByStringPatternString(shipCity: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByStringPatternString', { shipCity });
  }

  export async function postgresqlByBool(freightLess100: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByBool', { freightLess100 });
  }

  export async function postgresqlByDatetime(requiredDate: any): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByDatetime', { requiredDate });
  }

  export async function postgresqlByStringPatternDatetime(orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlByStringPatternDatetime', { orderDate });
  }

  export async function postgresqlOrders(employeeId: number, shipVia: string, freight: number, shipCountry: string, shipCity: string, freightLess1000: boolean, requiredDate: any, orderDate: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlOrders', { employeeId, shipVia, freight, shipCountry, shipCity, freightLess1000, requiredDate, orderDate });
  }

  export async function postgresqlProducts(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlProducts', {});
  }

  export async function postgresPerf(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresPerf', {});
  }

  export async function postgreScalarOutput(): Promise<number> {
    return await grok.data.query('@datagrok/dbtests:PostgreScalarOutput', {});
  }

  export async function postgresqlTableNormal(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlTableNormal', {});
  }

  export async function postgresqlTableWide(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlTableWide', {});
  }

  export async function postgresqlTableLong(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:PostgresqlTableLong', {});
  }

  export async function testForColumnsOnEmptyResult(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:TestForColumnsOnEmptyResult', {});
  }

  export async function simpleSelect(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SimpleSelect', {});
  }

  export async function snowflakeBinaryType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeBinaryType', {});
  }

  export async function snowflakeDateTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeDateTypes', {});
  }

  export async function snowflakeGeoType(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeGeoType', {});
  }

  export async function snowflakeNumericTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeNumericTypes', {});
  }

  export async function snowflakeSemiStructuredTypes(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeSemiStructuredTypes', {});
  }

  export async function snowflakePatternsAll(): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakePatternsAll', {});
  }

  export async function snowflakeIntTypePatternNone(id: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeIntTypePatternNone', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMore(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpMore', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMoreEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpMoreEq', { id });
  }

  export async function snowflakeStringTypeIntPatternOpLessEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpLessEq', { id });
  }

  export async function snowflakeStringTypeIntPatternOpLess(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpLess', { id });
  }

  export async function snowflakeStringTypeIntPatternOpIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpIn', { id });
  }

  export async function snowflakeStringTypeIntPatternOpNotIn(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpNotIn', { id });
  }

  export async function snowflakeStringTypeIntPatternOpMinMax(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpMinMax', { id });
  }

  export async function snowflakeStringTypeIntPatternOpNotEq(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypeIntPatternOpNotEq', { id });
  }

  export async function snowflakeDoubleTypePatternNone(some_number: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeDoubleTypePatternNone', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpMore(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternDoubleOpMore', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpMoreEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternDoubleOpMoreEq', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpLess(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternDoubleOpLess', { some_number });
  }

  export async function snowflakeStringTypePatternDoubleOpLessEq(some_number: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternDoubleOpLessEq', { some_number });
  }

  export async function snowflakeStringTypePatternStringOpContains(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternStringOpContains', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpStartsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternStringOpStartsWith', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpEndsWith(first_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternStringOpEndsWith', { first_name });
  }

  export async function snowflakeStringTypePatternStringOpIn(country: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternStringOpIn', { country });
  }

  export async function snowflakeStringTypePatternStringOpRegex(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakeStringTypePatternStringOpRegex', { email });
  }

  export async function snowflakePatternsAllParams(first_name: string, id: string, bool: boolean, email: string, some_number: string, country: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dbtests:SnowflakePatternsAllParams', { first_name, id, bool, email, some_number, country, date });
  }
}

export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/dbtests:Info', {});
  }

  export async function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): Promise<any> {
    return await grok.functions.call('@datagrok/dbtests:ExpectTable', { actual, expected });
  }
}
