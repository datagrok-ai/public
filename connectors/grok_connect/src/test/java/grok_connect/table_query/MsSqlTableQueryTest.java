package grok_connect.table_query;

import grok_connect.providers.MsSqlDataProvider;

public class MsSqlTableQueryTest extends TableQueryTest {
    public MsSqlTableQueryTest() {
        super(new MsSqlDataProvider());
    }

    @Override
    public String[] getExpectedEmptySchema() {
        return new String[] {"SELECT", "*", "FROM", "[events]"};
    }

    @Override
    public String[] getExpectedEmptySchemaTableWithDot() {
        return new String[] {"SELECT", "*", "FROM", "[public].[events]"};
    }

    @Override
    public String[] getExpectedEmptyFields() {
        return new String[] {"SELECT", "*", "FROM", "[public].[events]"};
    }

    @Override
    public String[] getExpectedEmptyFieldsLimit() {
        return new String[] {"SELECT", "top 50", "*", "FROM", "[public].[events]"};
    }

    @Override
    public String[] getExpectedFieldsWithoutDotLimit() {
        return new String[] {"SELECT", "top 50", "[id],", "[friendly_name],", "[name],", "[source],", "[session_id],",
                "[status],", "[description],", "[error_message],", "[error_stack_trace],",
                "[event_type_id]", "FROM", "[public].[events]"};
    }

    @Override
    public String[] getExpectedFieldsWithDotLimit() {
        return new String[] {"SELECT", "top 50", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[session_id],", "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],",
                "[events].[event_type_id]", "FROM", "[public].[events]"};
    }

    @Override
    public String[] getExpectedAggregationWithoutPattern() {
        return new String[] {"SELECT", "count(*) as [count(*)]", "FROM", "[public].[events]"};
    }

    @Override
    public String[] getExpectedAggregationWithPattern() {
        return new String[] {"SELECT", "sum([session_id]) as [sum(session_id)]", "FROM", "[public].[events]"};
    }

    @Override
    public String[] getAggregationAndGroupByLimitWithoutDot() {
        return new String[] {"SELECT", "top 50", "[id],", "[friendly_name],", "[name],", "[source],",
                "[status],", "[description],", "[error_message],", "[error_stack_trace],",
                "[event_type_id],", "count([session_id]) as [count(session_id)]", "FROM", "[public].[events]", "GROUP BY", "[id], [friendly_name], [name], [source], " +
                "[status], [description], [error_message], [error_stack_trace], [event_type_id]"};
    }

    @Override
    public String[] getAggregationAndGroupByLimitWithDot() {
        return new String[] {"SELECT", "top 50", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],", "count([events].[session_id]) as [count(events.session_id)]",
                "FROM", "[public].[events]", "GROUP BY", "[events].[id], [events].[friendly_name], [events].[name], [events].[source], " +
                "[events].[status], [events].[description], [events].[error_message], [events].[error_stack_trace], [events].[event_type_id]"};
    }

    @Override
    public String[] getAggregationAndGroupByAndHavingLimitWithoutDot() {
        return new String[] {"--input: string source", "SELECT", "top 50", "[id],", "[friendly_name],", "[name],", "[source],",
                "[status],", "[description],", "[error_message],", "[error_stack_trace],",
                "[event_type_id],", "count([session_id]) as [count(session_id)]", "FROM", "[public].[events]", "GROUP BY", "[id], [friendly_name], [name], [source], " +
                "[status], [description], [error_message], [error_stack_trace], [event_type_id]",
                "HAVING", "\t(@source([source]))"};
    }

    @Override
    public String[] getAggregationAndGroupByAndHavingLimitWithDot() {
        return new String[] {"--input: string events_source", "SELECT", "top 50", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],", "count([events].[session_id]) as [count(events.session_id)]",
                "FROM", "[public].[events]", "GROUP BY", "[events].[id], [events].[friendly_name], [events].[name], [events].[source], " +
                "[events].[status], [events].[description], [events].[error_message], [events].[error_stack_trace], [events].[event_type_id]",
                "HAVING", "\t(@events_source([events].[source]))"};
    }

    @Override
    public String[] getSimpleJoin() {
        return new String[] {"SELECT", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[session_id],", "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],",
                "[event_types].[id] as [event_types_id],", "[event_types].[friendly_name] as [event_types_friendly_name],",
                "[event_types].[name] as [event_types_name],", "[event_types].[description] as [event_types_description],",
                "[event_types].[error_message] as [event_types_error_message],", "[event_types].[error_stack_trace] as [event_types_error_stack_trace]",
                "FROM", "[public].[events]", "left join [event_types] on [events].[event_type_id] = [event_types].[id]"};
    }

    @Override
    public String[] getSelfJoin() {
        return new String[] {"SELECT", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[session_id],", "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],",
                "[events_1].[id] as [events_1_id],", "[events_1].[friendly_name] as [events_1_friendly_name],",
                "[events_1].[name] as [events_1_name],", "[events_1].[source] as [events_1_source],",
                "[events_1].[session_id] as [events_1_session_id],", "[events_1].[status] as [events_1_status],",
                "[events_1].[description] as [events_1_description],", "[events_1].[error_message] as [events_1_error_message],",
                "[events_1].[error_stack_trace] as [events_1_error_stack_trace],", "[events_1].[event_type_id] as [events_1_event_type_id]",
                "FROM", "[public].[events]", "left join [events] as [events_1] on [events].[event_type_id] = [events_1].[event_type_id]"};
    }

    @Override
    public String[] getSeveralJoins() {
        return new String[] {"SELECT", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[session_id],", "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],",
                "[event_types].[id] as [event_types_id],", "[event_types].[friendly_name] as [event_types_friendly_name],",
                "[event_types].[name] as [event_types_name],", "[event_types].[description] as [event_types_description],",
                "[event_types].[error_message] as [event_types_error_message],", "[event_types].[error_stack_trace] as [event_types_error_stack_trace],",
                "[sessions].[id] as [sessions_id],", "[sessions].[user_id],", "[sessions].[started],", "[sessions].[ended],", "[sessions].[token],",
                "[u].[id] as [u_id],", "[u].[email],", "[u].[first_name],", "[u].[last_name],", "[u].[status] as [u_status]",
                "FROM", "[public].[events]",
                "left join [event_types] on [events].[event_type_id] = [event_types].[id]",
                "right join [users_sessions] as [sessions] on [events].[session_id] = [sessions].[id]",
                "inner join [users] as [u] on [sessions].[user_id] = [u].[id]"
        };
    }

    @Override
    public String[] getSeveralOnJoin() {
        return new String[] {"SELECT", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[session_id],", "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],",
                "[t].[id] as [t_id],", "[t].[friendly_name] as [t_friendly_name],",
                "[t].[name] as [t_name],", "[t].[description] as [t_description],",
                "[t].[error_message] as [t_error_message],", "[t].[error_stack_trace] as [t_error_stack_trace]",
                "FROM", "[public].[events]", "left join [event_types] as [t] on [events].[event_type_id] = [t].[id]",
                " AND ", "[events].[name] = [t].[name]"
        };
    }
}
