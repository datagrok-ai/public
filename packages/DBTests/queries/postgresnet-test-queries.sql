--name: PostgresAll
--connection: PostgresTest
--meta.testExpectedRows: 830
select * from orders;
--end

--name: PostgresByInt
--connection: PostgresTest
--input: int orderid = 10330
--meta.testExpectedRows: 1
select * from orders where orderid = @orderid;
--end

--name: PostgresByStringPatternInt
--connection: PostgresTest
--input: string shipVia = '>2' {pattern: int}
--meta.testExpectedRows: 255
SELECT * FROM orders WHERE @shipVia(shipVia);
--end

--name: PostgresByDouble
--connection: PostgresTest
--input: double freight = 100.1
--meta.testExpectedRows: 187
SELECT * FROM orders WHERE freight >= @freight;
--end

--name: PostgresByStringChoices
--connection: PostgresTest
--input: string shipCountry = 'France' {choices: ["France", "Germany", "USA", "Finland"]}
--meta.testExpectedRows: 77
SELECT * FROM orders WHERE shipCountry = @shipCountry;
--end

--name: PostgresByStringPatternString
--connection: PostgresTest
--input: string shipCity = 'contains ran' {pattern: string}
--meta.testExpectedRows: 33
SELECT * FROM orders WHERE @shipCity(shipCity);
--end

--name: PostgresByBool
--connection: PostgresTest
--input: bool freightLess100 = true
--meta.testExpectedRows: 643
SELECT * FROM orders WHERE ((freight < 100) OR NOT @freightLess100);
--end

--name: PostgresByDatetime
--connection: PostgresTest
--input: datetime requiredDate
SELECT * FROM orders WHERE requiredDate >= @requiredDate;
--end

--name: PostgresByStringPatternDatetime
--connection: PostgresTest
--input: string orderDate {pattern: datetime}
SELECT * FROM orders WHERE @orderDate(orderDate);
--end

--name: PostgresOrders
--friendlyName: Orders
--connection: PostgresTest
--input: int employeeId = 5
--input: string shipVia = "= 3" {pattern: int}
--input: double freight = 10.0
--input: string shipCountry = "France" {choices: Query("SELECT DISTINCT shipCountry FROM Orders")}
--input: string shipCity = "starts with r" {pattern: string}
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

--name: PostgresProducts
--friendlyName: Products
--connection: PostgresTest
--input: int ProductID = 7
select * from Products
--end
