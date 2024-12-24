package grok_connect.table_query;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.PatternMatcher;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

@TestInstance(TestInstance.Lifecycle.PER_CLASS)
abstract class TableQueryTest {
    public static final String TEST_SCHEME = "public";
    public static final String TEST_TABLE = "events";
    public static final String[] TEST_TABLE_FIELDS = {
            "id", "friendly_name", "name", "source", "session_id",
            "status", "description", "error_message", "error_stack_trace",
            "event_type_id"
    };

    public static final String JOIN_TEST_TABLE_1 = "event_types";
    public static final String[] JOIN_TEST_TABLE_1_FIELDS = {
            "id", "friendly_name", "name", "description",
            "error_message", "error_stack_trace"
    };
    public static final String JOIN_TEST_TABLE_2 = "users_sessions";
    public static final String[] JOIN_TEST_TABLE_2_FIELDS = {
            "id", "user_id", "started", "ended", "token"
    };
    public static final String JOIN_TEST_TABLE_3 = "users";
    public static final String[] JOIN_TEST_TABLE_3_FIELDS = {
            "id", "email", "first_name", "last_name", "status"
    };

    private final JdbcDataProvider provider;
    private TableQuery testTableQuery;

    public TableQueryTest(JdbcDataProvider provider) {
        this.provider = provider;
        GrokConnect.providerManager = new ProviderManager();
    }

    @BeforeEach
    public void restore() {
        testTableQuery = new TableQuery();
        testTableQuery.tableName = TEST_TABLE;
        testTableQuery.schema = TEST_SCHEME;
        testTableQuery.connection = new DataConnection();
        testTableQuery.connection.dataSource = provider.descriptor.type;
        testTableQuery.params = new ArrayList<>();
    }

    @Test
    public void testEmptySchema() {
        testTableQuery.schema = null;
        testQuery(testTableQuery, getExpectedEmptySchema());
    }

    public abstract String[] getExpectedEmptySchema();

    @Test
    public void testEmptySchemaTableWithDot() {
        testTableQuery.schema = null;
        testTableQuery.tableName = TEST_SCHEME + "." + TEST_TABLE;
        testQuery(testTableQuery, getExpectedEmptySchemaTableWithDot());
    }

    public abstract String[] getExpectedEmptySchemaTableWithDot();

    @Test
    public void testEmptyFields() {
        testQuery(testTableQuery, getExpectedEmptyFields());
    }

    public abstract String[] getExpectedEmptyFields();

    @Test
    public void testEmptyFieldsLimit() {
        testTableQuery.limit = 50;
        testQuery(testTableQuery, getExpectedEmptyFieldsLimit());
    }

    public abstract String[] getExpectedEmptyFieldsLimit();

    @Test
    public void testFieldsWithoutDotLimit() {
        testTableQuery.limit = 50;
        testTableQuery.fields = Arrays.stream(TEST_TABLE_FIELDS).collect(Collectors.toList());
        testQuery(testTableQuery, getExpectedFieldsWithoutDotLimit());
    }

    public abstract String[] getExpectedFieldsWithoutDotLimit();

    @Test
    public void testFieldsWithDotLimit() {
        testTableQuery.limit = 50;
        testTableQuery.fields = Arrays.stream(TEST_TABLE_FIELDS)
                .map(f -> TEST_TABLE + "." + f)
                .collect(Collectors.toList());
        testQuery(testTableQuery, getExpectedFieldsWithDotLimit());
    }

    public abstract String[] getExpectedFieldsWithDotLimit();

    @Test
    public void testAggregationWithoutPattern() {
        testTableQuery.aggregations.add(new GroupAggregation(Stats.TOTAL_COUNT, "*", "count(*)"));
        testQuery(testTableQuery, getExpectedAggregationWithoutPattern());
    }

    public abstract String[] getExpectedAggregationWithoutPattern();

    @Test
    public void testAggregationWithPattern() {
        testTableQuery.aggregations.add(new GroupAggregation(Stats.SUM, "session_id", "sum(session_id)"));
        testQuery(testTableQuery, getExpectedAggregationWithPattern());
    }

    public abstract String[] getExpectedAggregationWithPattern();

    @Test
    public void testAggregationAndGroupByLimitWithoutDot() {
        testTableQuery.aggregations.add(new GroupAggregation(Stats.VALUE_COUNT, "session_id", "count(session_id)"));
        List<String> fields = Arrays.stream(TEST_TABLE_FIELDS).filter(f -> !f.equals("session_id")).collect(Collectors.toList());
        testTableQuery.fields = fields;
        testTableQuery.groupByFields = fields;
        testTableQuery.limit = 50;
        testQuery(testTableQuery, getAggregationAndGroupByLimitWithoutDot());
    }

    public abstract String[] getAggregationAndGroupByLimitWithoutDot();

    @Test
    public void testAggregationAndGroupByLimitWithDot() {
        testTableQuery.aggregations.add(new GroupAggregation(Stats.VALUE_COUNT, "events.session_id", "count(events.session_id)"));
        List<String> fields = Arrays.stream(TEST_TABLE_FIELDS).filter(f -> !f.equals("session_id"))
                .map(f -> TEST_TABLE + "." + f)
                .collect(Collectors.toList());
        testTableQuery.fields = fields;
        testTableQuery.groupByFields = fields;
        testTableQuery.limit = 50;
        testQuery(testTableQuery, getAggregationAndGroupByLimitWithDot());
    }

    public abstract String[] getAggregationAndGroupByLimitWithDot();

    @Test
    public void testAggregationAndGroupByAndHavingLimitWithoutDot() {
        testTableQuery.aggregations.add(new GroupAggregation(Stats.VALUE_COUNT, "session_id", "count(session_id)"));
        List<String> fields = Arrays.stream(TEST_TABLE_FIELDS).filter(f -> !f.equals("session_id")).collect(Collectors.toList());
        testTableQuery.fields = fields;
        testTableQuery.groupByFields = fields;
        testTableQuery.limit = 50;
        Map<String, Object> options = new HashMap<>();
        options.put("op", PatternMatcher.IN);
        FieldPredicate fieldPredicate = new FieldPredicate("source", "in", "string");
        List<String> values = new ArrayList<>();
        options.put("values", values);
        values.add("func");
        values.add("query");
        values.add("script");
        fieldPredicate.matcher = new PatternMatcher(options, "source");
        testTableQuery.having.add(fieldPredicate);
        testQuery(testTableQuery, getAggregationAndGroupByAndHavingLimitWithoutDot());
    }

    public abstract String[] getAggregationAndGroupByAndHavingLimitWithoutDot();


    @Test
    public void testAggregationAndGroupByAndHavingLimitWithDot() {
        testTableQuery.aggregations.add(new GroupAggregation(Stats.VALUE_COUNT, "events.session_id", "count(events.session_id)"));
        List<String> fields = Arrays.stream(TEST_TABLE_FIELDS).filter(f -> !f.equals("session_id"))
                .map(f -> TEST_TABLE + "." + f)
                .collect(Collectors.toList());
        testTableQuery.fields = fields;
        testTableQuery.groupByFields = fields;
        testTableQuery.limit = 50;
        Map<String, Object> options = new HashMap<>();
        options.put("op", PatternMatcher.IN);
        FieldPredicate fieldPredicate = new FieldPredicate("events.source", "in", "string");
        List<String> values = new ArrayList<>();
        options.put("values", values);
        values.add("func");
        values.add("query");
        values.add("script");
        fieldPredicate.matcher = new PatternMatcher(options, "events.source");
        testTableQuery.having.add(fieldPredicate);
        testQuery(testTableQuery, getAggregationAndGroupByAndHavingLimitWithDot());
    }

    public abstract String[] getAggregationAndGroupByAndHavingLimitWithDot();

    @Test
    public void testSimpleJoin() {
        testTableQuery.fields = Arrays.stream(TEST_TABLE_FIELDS).map(f -> TEST_TABLE + "." + f).collect(Collectors.toList());
        testTableQuery.fields.addAll(Arrays.stream(JOIN_TEST_TABLE_1_FIELDS).map(f -> JOIN_TEST_TABLE_1 + "." + f).collect(Collectors.toList()));
        TableJoin tableJoin = new TableJoin();
        tableJoin.joinType = "left";
        tableJoin.leftTableName = TEST_TABLE;
        tableJoin.rightTableName = JOIN_TEST_TABLE_1;
        tableJoin.leftTableKeys = new ArrayList<String>(){{ add("event_type_id"); }};
        tableJoin.rightTableKeys = new ArrayList<String>(){{ add("id"); }};
        testTableQuery.joins = new ArrayList<TableJoin>(){{ add(tableJoin); }};
        testQuery(testTableQuery, getSimpleJoin());
    }

    public abstract String[] getSimpleJoin();


    @Test
    public void testSelfJoin() {
        String selfJoinAlias = "events_1";
        testTableQuery.fields = Arrays.stream(TEST_TABLE_FIELDS).map(f -> TEST_TABLE + "." + f).collect(Collectors.toList());
        testTableQuery.fields.addAll(Arrays.stream(TEST_TABLE_FIELDS).map(f -> selfJoinAlias + "." + f).collect(Collectors.toList()));
        TableJoin tableJoin = new TableJoin();
        tableJoin.joinType = "left";
        tableJoin.leftTableName = TEST_TABLE;
        tableJoin.rightTableName = TEST_TABLE;
        tableJoin.rightTableAlias = selfJoinAlias;
        tableJoin.leftTableKeys = new ArrayList<String>(){{ add("event_type_id"); }};
        tableJoin.rightTableKeys = new ArrayList<String>(){{ add("event_type_id"); }};
        testTableQuery.joins = new ArrayList<TableJoin>(){{ add(tableJoin); }};
        testQuery(testTableQuery, getSelfJoin());
    }

    public abstract String[] getSelfJoin();

    @Test
    public void testSeveralJoins() {
        String sessionsAlias = "sessions";
        String usersAlias = "u";

        testTableQuery.fields = Arrays.stream(TEST_TABLE_FIELDS).map(f -> TEST_TABLE + "." + f).collect(Collectors.toList());
        testTableQuery.fields.addAll(Arrays.stream(JOIN_TEST_TABLE_1_FIELDS).map(f -> JOIN_TEST_TABLE_1 + "." + f).collect(Collectors.toList()));
        testTableQuery.fields.addAll(Arrays.stream(JOIN_TEST_TABLE_2_FIELDS).map(f -> sessionsAlias + "." + f).collect(Collectors.toList()));
        testTableQuery.fields.addAll(Arrays.stream(JOIN_TEST_TABLE_3_FIELDS).map(f -> usersAlias + "." + f).collect(Collectors.toList()));

        TableJoin tableJoin1 = new TableJoin();
        tableJoin1.joinType = "left";
        tableJoin1.leftTableName = TEST_TABLE;
        tableJoin1.rightTableName = JOIN_TEST_TABLE_1;
        tableJoin1.leftTableKeys = new ArrayList<String>(){{ add("event_type_id"); }};
        tableJoin1.rightTableKeys = new ArrayList<String>(){{ add("id"); }};
        testTableQuery.joins = new ArrayList<TableJoin>(){{ add(tableJoin1); }};

        TableJoin tableJoin2 = new TableJoin();
        tableJoin2.joinType = "right";
        tableJoin2.leftTableName = TEST_TABLE;
        tableJoin2.rightTableName = JOIN_TEST_TABLE_2;
        tableJoin2.rightTableAlias = sessionsAlias;
        tableJoin2.leftTableKeys = new ArrayList<String>(){{ add("session_id"); }};
        tableJoin2.rightTableKeys = new ArrayList<String>(){{ add("id"); }};
        testTableQuery.joins.add(tableJoin2);

        TableJoin tableJoin3 = new TableJoin();
        tableJoin3.joinType = "inner";
        tableJoin3.leftTableName = sessionsAlias;
        tableJoin3.rightTableName = JOIN_TEST_TABLE_3;
        tableJoin3.rightTableAlias = usersAlias;
        tableJoin3.leftTableKeys = new ArrayList<String>(){{ add("user_id"); }};
        tableJoin3.rightTableKeys = new ArrayList<String>(){{ add("id"); }};
        testTableQuery.joins.add(tableJoin3);

        testQuery(testTableQuery, getSeveralJoins());
    }

    public abstract String[] getSeveralJoins();

    @Test
    public void testSeveralOnJoin() {
        String rightTableAlias = "t";
        testTableQuery.fields = Arrays.stream(TEST_TABLE_FIELDS).map(f -> TEST_TABLE + "." + f).collect(Collectors.toList());
        testTableQuery.fields.addAll(Arrays.stream(JOIN_TEST_TABLE_1_FIELDS).map(f -> rightTableAlias + "." + f).collect(Collectors.toList()));
        TableJoin tableJoin = new TableJoin();
        tableJoin.joinType = "left";
        tableJoin.leftTableName = TEST_TABLE;
        tableJoin.rightTableName = JOIN_TEST_TABLE_1;
        tableJoin.rightTableAlias = rightTableAlias;
        tableJoin.leftTableKeys = new ArrayList<String>(){{
            add("event_type_id");
            add("name");
        }};
        tableJoin.rightTableKeys = new ArrayList<String>(){{
            add("id");
            add("name");
        }};
        testTableQuery.joins = new ArrayList<TableJoin>(){{ add(tableJoin); }};
        testQuery(testTableQuery, getSeveralOnJoin());
    }

    public abstract String[] getSeveralOnJoin();

    private void testQuery(TableQuery tableQuery, String[] expected) {
        String actual = provider.queryTableSql(tableQuery);
        String expectedQuery = getExpectedQuery(expected);
        Assertions.assertEquals(expectedQuery, actual);
    }

    private String getExpectedQuery(String[] expected) {
        return Arrays.stream(expected).collect(Collectors.joining(System.lineSeparator()));
    }
}
