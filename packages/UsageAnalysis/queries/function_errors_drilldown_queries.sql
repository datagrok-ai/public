--name: EventByErrorMessageAndFriendlyName
--meta.cache: true
--input: string errorMessage
--input: string friendlyName
--connection: System:Datagrok
select * from events e
where error_message = @errorMessage
and friendly_name = @friendlyName;
--end


--name: UpdateEventsIsErrorComment
--meta.cache: true
--input: string errorMessage
--input: string friendlyName
--input: bool isError
--input: string comment
--connection: System:Datagrok
update events
set is_error = @isError,
comment = @comment
where error_message = @errorMessage
and friendly_name = @friendlyName;
--end