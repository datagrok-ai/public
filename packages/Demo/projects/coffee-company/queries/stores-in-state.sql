--name: stores in @state 
--input: string state = 'NY' { pattern: string }
--tags: unit-test
--meta.testExpected: 645
select * from starbucks_us where @state(state)
--end
