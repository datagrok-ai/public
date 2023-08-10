--name: stores in @state with sales greater than @sales
--input: string state = 'GA'
--input: int sales = '99581.29688'

select * from starbucks_us where state = @state and sales > @sales
--end
