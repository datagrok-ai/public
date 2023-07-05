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
SELECT pg_sleep(3);
--batch
SELECT max(float_data) FROM Test_Long;
--end

--name: PostgresqlScalarCacheTestInt
--friendlyName: PostgresqlScalarCacheTestInt
--connection: PostgreSQLDBTests
--output: int count
--meta.cache: true
SELECT pg_sleep(3);
--batch
SELECT count(*) FROM Test_Long;
--end

--name: PostgresqlScalarCacheTestString
--friendlyName: PostgresqlScalarCacheTestString
--connection: PostgreSQLDBTests
--output: string first_name
--meta.cache: true
SELECT pg_sleep(3);
--batch
SELECT first_name FROM MOCK_DATA WHERE id = 10;
--end

--name: PostgresqlCachedConnTest
--friendlyName: PostgresqlCachedConnTest
--connection: PostgreSQLDBTestsCached
SELECT pg_sleep(3);
--batch
SELECT * FROM MOCK_DATA;
--end

--name: PostgresqlCacheInvalidateQueryTest
--friendlyName: PostgresqlCacheInvalidateQueryTest
--connection: PostgreSQLDBTests
--meta.cache: true
--meta.invalidate: * * * * *
SELECT pg_sleep(3);
--batch
SELECT * FROM MOCK_DATA;
--end