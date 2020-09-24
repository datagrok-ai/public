
--name: #{NAME}
--connection: #{CONNECTION}
--input: int id
--output: dataframe result
select * from country where id = @id
--end
