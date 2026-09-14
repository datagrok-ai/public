--name: compounds for @target
--friendlyName: Browse | Compounds For Target
--description: Retrieves compounds.
--connection: Demo
--input: string target = "CHEMBL1827" [ChEMBL target identifier]
--test: Dbtests:expectTable(CompoundsForTarget(), OpenFile('x.d42'))
--meta.testExpectedRows: 5
--meta.cache: server
select 1
--end

--name: users
--connection: System:Datagrok
--output: dataframe result
select * from users
--end

--name: no connection
--input: int n
select 1
