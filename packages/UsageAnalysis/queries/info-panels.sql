--name: userById
--tags: panel
--render: RowToTable
--input: string id { semType: user_id }
--connection: datagrok
select u.first_name, u.last_name from users u
where id = @id
--end

