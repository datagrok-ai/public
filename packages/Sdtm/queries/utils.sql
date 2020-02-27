--name: Studies
--input: string race {pattern: string}
--connection: clin
select distinct studyid as study
from (
  select lb.studyid as studyid
  from lb
  inner join dm on lb.studyid = dm.studyid and lb.usubjid = dm.usubjid
  where @race(race)
) as t
order by study
--end


--name: Tests
--input: string race {pattern: string}
--input: string studies {pattern: string}
--connection: clin
select distinct lbtest as test
from (
  select lbtest
  from lb
  inner join dm on lb.studyid = dm.studyid and lb.usubjid = dm.usubjid
  where @studies(lb.studyid) and @race(race)
) as t
order by test
--end


--name: Methods
--input: string race {pattern: string}
--input: string studies {pattern: string}
--input: string tests {pattern: string}
--connection: clin
select distinct lbmethod as method
from (
  select lbmethod
  from lb
  inner join dm on lb.studyid = dm.studyid and lb.usubjid = dm.usubjid
  where @race(race) and @studies(lb.studyid) and @tests(lbtest)
) as t
where lbmethod != ''
order by method
--end
