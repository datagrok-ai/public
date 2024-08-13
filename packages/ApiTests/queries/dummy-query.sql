--name: dummyPackageQuery
--connection: System:Datagrok
--input: double x
SELECT id, description from event_types;
--end
