--name: stores in @state with sales greater than @sales
--input: string state = 'NY' 
--input: int sales
select * from starbucks_us where state = @state and sales > @sales
--end
