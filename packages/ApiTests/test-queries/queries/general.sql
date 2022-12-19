--name: WithOption
--connection: PostgresNorthwind
--meta.cache: false
select 1;
--end

--name: BatchMode
--connection: PostgresNorthwind
--meta.batchMode: true
select 1;
--batch
select 2;
--end
