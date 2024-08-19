--input: string when {pattern: datetime}
--meta.search: new users @when
--connection: System:Datagrok
select * from users
where @when(joined)
