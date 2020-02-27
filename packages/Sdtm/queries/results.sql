--name: Results
--input: string studies {pattern: string}
--input: string tests {pattern: string}
--input: string sex {pattern: string}
--input: string race {pattern: string}
--input: string method {pattern: string}
--input: string limit {pattern: int}
--connection: clin
select
  lb.studyid as studyid,
  lb.usubjid as usubjid,
  lbdy,
  lbtest as test,
  lbstresn,
  lbstresu as units,
  lbmethod as method,
  sex,
  race
from lb
inner join dm on lb.studyid = dm.studyid and lb.usubjid = dm.usubjid
where
  @studies(lb.studyid) and @tests(lbtest) and @method(lbmethod) and @sex(sex) and @race(race)
group by
  lb.studyid, lb.usubjid, lbdy, lbtest, lbstresn, lbmethod, lbstresu, sex, race
limit @limit
--end


--name: ResultsCount
--input: string studies {pattern: string}
--input: string tests {pattern: string}
--input: string sex {pattern: string}
--input: string race {pattern: string}
--input: string method {pattern: string}
--connection: clin
select count(*) from (
  select
    lb.studyid as studyid,
    lb.usubjid as usubjid,
    lbdy,
    lbtest as test,
    lbstresn,
    lbstresu as units,
    lbmethod as method,
    sex,
    race
  from lb
  inner join dm on lb.studyid = dm.studyid and lb.usubjid = dm.usubjid
  where
    @studies(lb.studyid) and @tests(lbtest) and @method(lbmethod) and @sex(sex) and @race(race)
  group by
    lb.studyid, lb.usubjid, lbdy, lbtest, lbstresn, lbmethod, lbstresu, sex, race
) as t
--end
