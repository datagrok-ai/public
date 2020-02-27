--name: MedicalHistory
--input: string study
--input: string subj
--connection: clin
select mhterm, mhenrf
from mh
where studyid = @study and usubjid = @subj and mhoccur = 'Y'
order by mhterm
--end
