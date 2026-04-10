--name: DatabricksBenchSmall
--friendlyName: DatabricksBenchSmall
--connection: DatabricksDBTests
SELECT * FROM workspace.grok_bench.bench_small
--end

--name: DatabricksBenchMedium
--friendlyName: DatabricksBenchMedium
--connection: DatabricksDBTests
SELECT * FROM workspace.grok_bench.bench_medium
--end

--name: DatabricksBenchLarge
--friendlyName: DatabricksBenchLarge
--connection: DatabricksDBTests
SELECT * FROM workspace.grok_bench.bench_large
--end

--name: DatabricksBenchCategorical
--friendlyName: DatabricksBenchCategorical
--connection: DatabricksDBTests
SELECT * FROM workspace.grok_bench.bench_categorical
--end

--name: DatabricksSelect1
--friendlyName: DatabricksSelect1
--connection: DatabricksDBTests
SELECT 1
--end
