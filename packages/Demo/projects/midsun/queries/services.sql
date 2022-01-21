select
  s.date, s.durationMin, s.pathDurationMin,
  st.name as servicetype,
  e.name as engineer
from service s
  join servicetype st on s.ServiceTypeId = st.id
  join engineertoservice ets on s.id = ets.ServiceId
  join engineer e on ets.EngineerId = e.id