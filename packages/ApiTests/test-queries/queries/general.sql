--name: WithOption
--connection: PostgresTest
--meta.cache: false
select 1;
--end

--name: BatchMode
--connection: PostgresTest
--meta.batchMode: true
select 1;
--batch
select 2;
--end
