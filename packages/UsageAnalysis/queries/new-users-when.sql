--input: string when {pattern: datetime}
--meta.search: new users @when
select * from users
where @when(joined)