--name: PostgresqlTestCacheTableWide
--friendlyName: PostgresqlTestCacheTableWide
--connection: PostgreSQLApiTests
--meta.cache: true
SELECT * FROM Test_Wide;
--end

--name: PostgresqlTestCacheTableNormal
--friendlyName: PostgresqlTestCacheTableNormal
--connection: PostgreSQLApiTests
--meta.cache: true
SELECT * FROM Test_Normal;
--end
