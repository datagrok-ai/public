--name: userById
--tags: panel
--render: RowToTable
--input: string id { semType: user_id }
--connection: System:DatagrokAdmin
select u.first_name, u.last_name from users u
where id = @id
--end

--name: userInfoByEmailPanel
--connection: System:DatagrokAdmin
--tags: panel
--input: string email {semType:email}

select u.first_name, u.last_name, u.joined, u.email
from users u
where email = @email
--end