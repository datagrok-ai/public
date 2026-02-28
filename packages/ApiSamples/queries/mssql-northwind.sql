--name: MSSQLByInt
--connection: MSSQLNorthwind
--input: int orderid = 10
select * from orders where orderid = @orderid;
--end
