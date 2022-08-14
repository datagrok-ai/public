--name: WithOption
--connection: PostgresNetNorthwind
--meta.cache: false
select 1;
--end

--name: BatchMode
--connection: PostgresNetNorthwind
--meta.batchMode: true
select 1;
--batch
select 2;
--end
