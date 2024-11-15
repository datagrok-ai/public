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
        return new String[] {"SELECT", "top 50", "[id],", "[friendly_name],", "[name],", "[source],",
                "[status],", "[description],", "[error_message],", "[error_stack_trace],",
                "[event_type_id],", "count([session_id]) as [count(session_id)]", "FROM", "[public].[events]", "GROUP BY", "[id], [friendly_name], [name], [source], " +
                "[status], [description], [error_message], [error_stack_trace], [event_type_id]",
                "HAVING", "\t((LOWER([source]) IN ('func','query','script')))"};
    }

    @Override
    public String[] getAggregationAndGroupByAndHavingLimitWithDot() {
        return new String[] {"SELECT", "top 50", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],", "count([events].[session_id]) as [count(events.session_id)]",
                "FROM", "[public].[events]", "GROUP BY", "[events].[id], [events].[friendly_name], [events].[name], [events].[source], " +
                "[events].[status], [events].[description], [events].[error_message], [events].[error_stack_trace], [events].[event_type_id]",
                "HAVING", "\t((LOWER([events].[source]) IN ('func','query','script')))"};
    }

    @Override
    public String[] getSimpleJoin() {
        return new String[] {"SELECT", "[events].[id],", "[events].[friendly_name],", "[events].[name],", "[events].[source],",
                "[events].[session_id],", "[events].[status],", "[events].[description],", "[events].[error_message],",
                "[events].[error_stack_trace],", "[events].[event_type_id],",
                "[event_types].[id] as [event_types_id],", "[event_types].[friendly_name] as [event_types_friendly_name],",
                "[event_types].[name] as [event_types_name],", "[event_types].[description] as [event_types_description],",
                "[event_types].[error_message] as [event_types_error_message],", "[event_types].[error_stack_trace] as [event_types_error_stack_trace]",
                "FROM", "[public].[events]", "left join event_types on [events].[event_type_id] = [event_types].[id]"};
    }
}
