--name: Tests
--input: string studies {pattern: string}
select distinct
  LBTEST as Test
from
  LB
where
  @studies(STUDYID)
--end
