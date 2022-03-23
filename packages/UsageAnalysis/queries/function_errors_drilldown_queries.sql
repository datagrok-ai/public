--name: EventByErrorMessageAndFriendlyName
--input: string errorMessage
--input: string friendlyName
--connection: System:DatagrokAdmin
select * from events e
where error_message = @errorMessage
and friendly_name = @friendlyName;
--end


--name: UpdateEventsIsErrorComment
--input: string errorMessage
--input: string friendlyName
--input: bool isError
--input: string comment
--connection: System:DatagrokAdmin
update events
set is_error = @isError,
comment = @comment
where error_message = @errorMessage
and friendly_name = @friendlyName;
--end