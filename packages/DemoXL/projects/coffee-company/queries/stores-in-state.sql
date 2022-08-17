--name: stores in @state 
--input: string state = 'NY' { pattern: string }
--meta.testExpected: 645
select * from starbucks_us where @state(state)
--end
