--name: dummyPackageQuery
--connection: System:Datagrok
--input: double x
SELECT id, description, @x as res from event_types limit 10;
--end
