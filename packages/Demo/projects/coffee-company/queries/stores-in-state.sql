--name: stores in @state
--input: string state { pattern: string }
select * from starbucks_us where @state(state)
--end
