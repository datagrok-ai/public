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
