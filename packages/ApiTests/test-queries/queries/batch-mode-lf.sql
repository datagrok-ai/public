--name: BatchModeLf
--connection: PostgresTest
--meta.batchMode: true
select 1;
--batch
select 2;
--end
