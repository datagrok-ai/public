--name: PostgresqlAll
--connection: PostgreSQLTest
--meta.testExpectedRows: 830
select * from orders;
--end

--name: PostgresqlByInt
--connection: PostgreSQLTest
--input: int orderid = 10330
--meta.testExpectedRows: 1
select * from orders where orderid = @orderid;
--end

--name: PostgresqlByStringPatternInt
--connection: PostgreSQLTest
--input: string shipVia = '>2' {pattern: int}
--meta.testExpectedRows: 255
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: PostgresqlByDouble
--connection: PostgreSQLTest
--input: double freight = 100.1
--meta.testExpectedRows: 187
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: PostgresqlByStringChoices
--connection: PostgreSQLTest
--input: string shipCountry = 'France' {choices: ["France", "Germany", "USA", "Finland"]}
--meta.testExpectedRows: 77
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: PostgresqlByStringPatternString
--connection: PostgreSQLTest
--input: string shipCity = 'contains ran' {pattern: string}
--meta.testExpectedRows: 33
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: PostgresqlByBool
--connection: PostgreSQLTest
--input: bool freightLess100 = true
--meta.testExpectedRows: 643
SELECT * FROM orders WHERE ((freight < 100) OR NOT @freightLess100);
--end

--name: PostgresqlByDatetime
--connection: PostgreSQLTest
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: PostgresqlByStringPatternDatetime
--connection: PostgreSQLTest
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end

--name: PostgresqlOrders
--friendlyName: Orders
--connection: PostgreSQLTest
--tags: unit-test
--input: int employeeId = 5
--input: string shipVia = "3" {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCity = "contains r" {pattern: string}
--input: bool freightLess1000 = true
--input: datetime requiredDate = "1/1/1995"
--input: string orderDate = "after 1/1/1995" {pattern: datetime}

SELECT * FROM Orders WHERE (employeeId = @employeeId)
    AND (freight >= @freight)
    AND @shipVia(shipVia)
    AND ((freight < 1000) OR NOT @freightLess1000)
    AND (shipCountry = @shipCountry)
    AND @shipCity(shipCity)
    AND @orderDate(orderDate)
    AND (requiredDate >= @requiredDate)
--end

--name: PostgresqlProducts
--friendlyName: Products
--connection: PostgreSQLTest
--meta.testExpectedRows: 77
select * from Products
--end

--name: PostgresPerf
--friendlyName: Perf
--connection: PostgresTest
--output: dataframe res
SELECT * FROM generate_series(1,200000);
--end

--name: PostgreScalarOutput
--connection: PostgreSQLTest
--output: int result
select Count(*) from orders;
--end

--friendlyName: ChemblPerfGenerated
--connection: PostgresChemblTest
--input: int num = 5000000
select
    left(md5(i::text), 10),
    md5(random()::text),
    md5(random()::text),
    md5(random()::text),
    left(md5(random()::text), 4)
from generate_series(1, @num) s(i)
--end

--friendlyName: Compounds
--connection: PostgresChemblTest
select * from compound_structures
--end

--name: PostgresqlTableNormal
--friendlyName: PostgresqTableNormal
--connection: PostgreSQLDBTests
SELECT * FROM Test_Normal;
--end

--name: PostgresqlTableWide
--friendlyName: PostgresqlTableWide
--connection: PostgreSQLDBTests
SELECT * FROM Test_Wide;
--end

--name: PostgresqlTableLong
--friendlyName: PostgresqlTableLong
--connection: PostgreSQLDBTests
SELECT * FROM Test_Long;
--end

--name: TestForColumnsOnEmptyResult
--friendlyName: TestForColumnsOnEmptyResult
--connection: PostgreSQLDBTests
SELECT * FROM mock_data WHERE id > 30;
--end

--name: SimpleSelect
--friendlyName: SimpleSelect
--connection: PostgreSQLDBTests
SELECT 1;
--end