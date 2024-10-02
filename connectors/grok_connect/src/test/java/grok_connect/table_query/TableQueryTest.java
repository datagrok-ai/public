package grok_connect.table_query;

import grok_connect.connectors_info.DataConnection;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.PatternMatcher;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
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
class TableQueryTest {
    private static final String TEST_TABLE = "users";
    private static final String TEST_SCHEME = "public";
    private JdbcDataProvider postgres;
    private JdbcDataProvider mssql;
    private TableQuery limitTest;
    private TableQuery aggregationTest;

    @BeforeAll
    public void init() {
        ProviderManager providerManager = new ProviderManager();
        postgres = providerManager.getByName("Postgres");
        mssql = providerManager.getByName("MS SQL");
    }

    @BeforeEach
    public void restore() {
        limitTest = new TableQuery();
        limitTest.tableName = TEST_TABLE;
        limitTest.schema = TEST_SCHEME;
        List<String> fields = new ArrayList<>();
        fields.add("first_name");
        fields.add("last_name");
        fields.add("email");
        fields.add("address");
        limitTest.fields = fields;
        limitTest.limit = 50;
        limitTest.connection = new DataConnection();
        aggregationTest = new TableQuery();
        aggregationTest.connection = new DataConnection();
        aggregationTest.tableName = TEST_TABLE;
        aggregationTest.schema = TEST_SCHEME;
    }

    @Test
    public void testLimitAtEnd() {
        limitTest.connection.dataSource = postgres.descriptor.type;
        String sqlQuery = postgres.queryTableSql(limitTest.connection, limitTest);
        String[] expected = new String[] {"SELECT", "\"first_name\",", "\"last_name\",", "\"email\",", "\"address\"", "FROM",
                "\"public\".\"users\"", "limit 50"};
        Assertions.assertEquals(Arrays.stream(expected).collect(Collectors.joining(System.lineSeparator())),
                sqlQuery);
    }

    @Test
    public void testLimitAtStart() {
        limitTest.connection.dataSource = mssql.descriptor.type;
        String sqlQuery = mssql.queryTableSql(limitTest.connection, limitTest);
        String[] expected = new String[] {"SELECT", "top 50", "[first_name],", "[last_name],", "[email],", "[address]", "FROM",
                "[public].[users]"};
        Assertions.assertEquals(getExpectedQuery(expected),
                sqlQuery);
    }

    @Test
    public void testAggregation() {
        aggregationTest.connection.dataSource = postgres.descriptor.type;
        List<GroupAggregation> aggregations = new ArrayList<>();
        aggregations.add(new GroupAggregation(Stats.TOTAL_COUNT, "*", "count(*)"));
        aggregationTest.aggregations = aggregations;
        DataConnection connection = new DataConnection();
        connection.dataSource = postgres.descriptor.type;
        aggregationTest.connection = connection;
        String sqlQuery = postgres.queryTableSql(connection, aggregationTest);
        String[] expected = new String[] {"SELECT", "count(*) as \"count(*)\"", "FROM", "\"public\".\"users\""};
        Assertions.assertEquals(getExpectedQuery(expected), sqlQuery);
    }

    @Test
    public void testAggregationAndGroupBy() {
        aggregationTest.connection.dataSource = postgres.descriptor.type;
        List<GroupAggregation> aggregations = new ArrayList<>();
        aggregations.add(new GroupAggregation(Stats.TOTAL_COUNT, "*", "count(*)"));
        aggregationTest.aggregations = aggregations;
        List<String> fields = new ArrayList<>();
        fields.add("country");
        fields.add("region");
        aggregationTest.fields = fields;
        List<String> groupByFields = new ArrayList<>();
        groupByFields.add("country");
        groupByFields.add("region");
        aggregationTest.groupByFields = groupByFields;
        String sqlQuery = postgres.queryTableSql(aggregationTest.connection, aggregationTest);
        String[] expected = new String[] {"SELECT",  "\"country\",", "\"region\",", "count(*) as \"count(*)\"", "FROM", "\"public\".\"users\"",
                "GROUP BY", "\"country\", \"region\""};
        Assertions.assertEquals(getExpectedQuery(expected), sqlQuery);
    }

    @Test
    public void testAggregationAndGroupByAndHavingLimitAtEnd() {
        prepareComplexQuery(postgres.descriptor.type);
        String sqlQuery = postgres.queryTableSql(aggregationTest.connection, aggregationTest);
        String[] expected = new String[] {"SELECT", "\"country\",", "\"region\",", "count(*) as \"count(*)\"", "FROM", "\"public\".\"users\"",
                "GROUP BY", "\"country\", \"region\"", "HAVING", "\t((LOWER(country) IN ('spain','ukraine','brazil')))", "limit 50"};
        Assertions.assertEquals(getExpectedQuery(expected), sqlQuery);
    }

    @Test
    public void testAggregationAndGroupByAndHavingLimitAtStart() {
        prepareComplexQuery(mssql.descriptor.type);
        String sqlQuery = mssql.queryTableSql(aggregationTest.connection, aggregationTest);
        String[] expected = new String[] {"SELECT", "top 50", "[country],", "[region],", "count(*) as [count(*)]", "FROM", "[public].[users]",
                "GROUP BY", "[country], [region]", "HAVING", "\t((LOWER(country) IN ('spain','ukraine','brazil')))"};
        Assertions.assertEquals(getExpectedQuery(expected), sqlQuery);
    }

    private void prepareComplexQuery(String providerType) {
        aggregationTest.connection.dataSource = providerType;
        List<GroupAggregation> aggregations = new ArrayList<>();
        aggregations.add(new GroupAggregation(Stats.TOTAL_COUNT, "*", "count(*)"));
        aggregationTest.aggregations = aggregations;
        List<String> fields = new ArrayList<>();
        fields.add("country");
        fields.add("region");
        aggregationTest.fields = fields;
        List<String> groupByFields = new ArrayList<>();
        groupByFields.add("country");
        groupByFields.add("region");
        aggregationTest.groupByFields = groupByFields;
        FieldPredicate fieldPredicate = new FieldPredicate("country", "in", "string");
        Map<String, Object> options = new HashMap<>();
        options.put("op", PatternMatcher.IN);
        List<String> values = new ArrayList<>();
        options.put("values", values);
        values.add("Spain");
        values.add("Ukraine");
        values.add("Brazil");
        fieldPredicate.matcher = new PatternMatcher(options, "country");
        List<FieldPredicate> having = new ArrayList<>();
        having.add(fieldPredicate);
        aggregationTest.having = having;
        aggregationTest.limit = 50;
    }

    private String getExpectedQuery(String[] expected) {
        return Arrays.stream(expected).collect(Collectors.joining(System.lineSeparator()));
    }
}
