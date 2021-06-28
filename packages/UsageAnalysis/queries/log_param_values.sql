--name: Log Parameter Values
--connection: System:Datagrok

select et.name as event_name, e.event_time, ep.name as param_name, epv.value 
from event_types et, events e, event_parameters ep, event_parameter_values epv
where e.event_type_id = et.id 
  and epv.event_id = e.id
  and epv.parameter_id = ep.id
  and et.name not like 'error%'