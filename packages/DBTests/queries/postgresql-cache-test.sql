--name: PostgresqlTestCacheTableWide
--friendlyName: PostgresqlTestCacheTableWide
--connection: PostgreSQLDBTests
--meta.cache: true
SELECT * FROM Test_Wide;
--end

--name: PostgresqlTestCacheTableNormal
--friendlyName: PostgresqlTestCacheTableNormal
--connection: PostgreSQLDBTests
--meta.cache: true
SELECT * FROM Test_Normal;
--end

--name: PostgresqlScalarCacheTestFloat
--friendlyName: PostgresqlScalarCacheTestFloat
--connection: PostgreSQLDBTests
--output: double maxValue
--meta.cache: true
SELECT max(float_data) as maxValue, pg_sleep(3) FROM Test_Long;
--end

--name: PostgresqlScalarCacheTestInt
--friendlyName: PostgresqlScalarCacheTestInt
--connection: PostgreSQLDBTests
--output: int count
--meta.cache: true
SELECT count(*) as count, pg_sleep(3) FROM Test_Long;
--end

--name: PostgresqlScalarCacheTestString
--friendlyName: PostgresqlScalarCacheTestString
--connection: PostgreSQLDBTests
--output: string first_name
--meta.cache: true
SELECT first_name, pg_sleep(3) FROM MOCK_DATA WHERE id = 10;
--end

--name: PostgresqlScalarCacheTestDate
--friendlyName: PostgresqlScalarCacheTestDate
--connection: PostgreSQLDBTests
--output: datetime date
--meta.cache: true
SELECT date, pg_sleep(3) FROM MOCK_DATA WHERE id = 10;
--end

--name: PostgresqlCachedConnTest
--friendlyName: PostgresqlCachedConnTest
--connection: PostgreSQLDBTestsCached
SELECT *, pg_sleep(0.1) FROM MOCK_DATA;
--end

--name: PostgresqlCacheInvalidateQueryTest
--friendlyName: PostgresqlCacheInvalidateQueryTest
--connection: PostgreSQLDBTests
--meta.cache: true
SELECT *, pg_sleep(0.1) FROM MOCK_DATA;
--end

--name: PostgresqlScalarCacheInvalidationTestFloat
--friendlyName: PostgresqlScalarCacheInvalidationTestFloat
--connection: PostgreSQLDBTests
--output: double maxValue
--meta.cache: true
--meta.invalidate: 0 1 * * *
SELECT max(float_data) as maxValue, pg_sleep(3) FROM Test_Long;
--end

--name: PostgresqlScalarCacheInvalidationTestInt
--friendlyName: PostgresqlScalarCacheInvalidationTestInt
--connection: PostgreSQLDBTests
--output: int count
--meta.cache: true
--meta.invalidate: 0 1 * * *
SELECT count(*) as count, pg_sleep(3) FROM Test_Long;
--end

--name: PostgresqlScalarCacheInvalidationTestString
--friendlyName: PostgresqlScalarCacheInvalidationTestString
--connection: PostgreSQLDBTests
--output: string first_name
--meta.cache: true
--meta.invalidate: 0 1 * * *
SELECT first_name FROM MOCK_DATA, pg_sleep(3) WHERE id = 10;
--end

--name: PostgresqlScalarCacheInvalidationTestDate
--friendlyName: PostgresqlScalarCacheInvalidationTestDate
--connection: PostgreSQLDBTests
--output: datetime date
--meta.cache: true
--meta.invalidate: 0 1 * * *
SELECT date FROM MOCK_DATA, pg_sleep(3) WHERE id = 10;
--end
