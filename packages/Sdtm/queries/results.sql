--name: ResultsPreview
--input: string studies {pattern: string}
--input: string tests {pattern: string}
--input: int limit = 100
select
  STUDYID,
  USUBJID,
  LBDY,
  LBTEST,
  LBORRES
from
  LB
where
  @studies(STUDYID) and @tests(LBTEST)
group by
  STUDYID, USUBJID, LBDY, LBTEST, LBORRES
limit @limit
--end


--name: Results
--input: string studies {pattern: string}
--input: string tests {pattern: string}
select
  STUDYID,
  USUBJID,
  LBDY,
  LBTEST,
  LBORRES
from
  LB
where
  @studies(STUDYID) and @tests(LBTEST)
group by
  STUDYID, USUBJID, LBDY, LBTEST, LBORRES
--end
