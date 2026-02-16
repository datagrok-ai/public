package grok_connect.utils;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import serialization.DataFrame;
import serialization.StringColumn;

import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

public class SqlAnnotatorTest {

    @BeforeAll
    static void setup() {
        GrokConnect.providerManager = new ProviderManager();
    }

    private DataQuery query(String sql, String dataSource) {
        DataQuery q = new DataQuery();
        q.query = sql;

        q.connection = new DataConnection();
        q.connection.id = "conn-1";
        q.connection.dataSource = dataSource;

        return q;
    }

    private DataFrame df(String... columnNames) {
        DataFrame df = new DataFrame();
        for (String c : columnNames)
            df.addColumn(new StringColumn(c));
        return df;
    }

    static class SqlCase {
        final String sql;
        final String dataSource;
        final String[] columns;
        final ColumnAssertion assertion;

        SqlCase(String sql, String dataSource, String[] columns, ColumnAssertion assertion) {
            this.sql = sql;
            this.dataSource = dataSource;
            this.columns = columns;
            this.assertion = assertion;
        }
    }

    interface ColumnAssertion {
        void assertColumn(serialization.Column<?> c);
    }

    static Stream<SqlCase> starQueries() {
        return Stream.of(

                new SqlCase(
                        "select * from schema1.table1",
                        "Postgres",
                        new String[]{"a", "b"},
                        c -> {
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertEquals(c.getName(), c.getTag(Tags.DbColumn));
                            assertEquals("schema1", c.getTag(Tags.DbSchema));
                        }
                )
        );
    }

    static Stream<SqlCase> explicitQueries() {
        return Stream.of(

                new SqlCase(
                        "select t.col1 from table1 t",
                        "Postgres",
                        new String[]{"col1"},
                        c -> {
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertEquals("col1", c.getTag(Tags.DbColumn));
                        }
                )
        );
    }

    static Stream<SqlCase> derivedQueries() {
        return Stream.of(

                new SqlCase(
                        "select col1 + col2 as sum from table1",
                        "Postgres",
                        new String[]{"sum"},
                        c -> {
                            assertEquals("true", c.getTag(Tags.IsDerived));
                            assertEquals("col1 + col2", c.getTag(Tags.DbExpression));
                        }
                )
        );
    }

    static Stream<SqlCase> joinQueries() {
        return Stream.of(

                new SqlCase(
                        "select a.id, b.name from table1 a join table2 b on a.id = b.id",
                        "Postgres",
                        new String[]{"id", "name"},
                        c -> {
                            assertNotNull(c.getTag(Tags.DbTable));
                            assertNotNull(c.getTag(Tags.DbColumn));
                        }
                ),

                new SqlCase(
                        "select a.id, b.name " +
                                "from table1 a left join table2 b on a.id = b.id",
                        "Postgres",
                        new String[]{"id", "name"},
                        c -> {
                            assertNotNull(c.getTag(Tags.DbTable));
                            assertNotNull(c.getTag(Tags.DbColumn));
                        }
                ),

                new SqlCase(
                        "select b.name " +
                                "from table1 a left join table2 b on a.id = b.id",
                        "Postgres",
                        new String[]{"name"},
                        c -> {
                            assertEquals("table2", c.getTag(Tags.DbTable));
                            assertEquals("name", c.getTag(Tags.DbColumn));
                        }
                )
        );
    }

    static Stream<SqlCase> quotedIdentifierQueries() {
        return Stream.of(

                new SqlCase(
                        "select \"Col1\" from \"MySchema\".\"MyTable\"",
                        "Postgres",
                        new String[]{"Col1"},
                        c -> {
                            assertEquals("MyTable", c.getTag(Tags.DbTable));
                            assertEquals("Col1", c.getTag(Tags.DbColumn));
                            assertEquals("MySchema", c.getTag(Tags.DbSchema));
                        }
                ),

                new SqlCase(
                        "select `col1` from `table1`",
                        "Databricks",
                        new String[]{"col1"},
                        c -> {
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertEquals("col1", c.getTag(Tags.DbColumn));
                        }
                )
        );
    }

    static Stream<SqlCase> mssqlQueries() {
        return Stream.of(

                new SqlCase(
                        "select [c1] from [dbo].[MyTable]",
                        "MS SQL",
                        new String[]{"c1"},
                        c -> {
                            assertEquals("MyTable", c.getTag(Tags.DbTable));
                            assertEquals("c1", c.getTag(Tags.DbColumn));
                            assertEquals("dbo", c.getTag(Tags.DbSchema));
                        }
                ),

                new SqlCase(
                        "select t.* from [dbo].[MyTable] t",
                        "MS SQL",
                        new String[]{"a", "b"},
                        c -> {
                            assertEquals("MyTable", c.getTag(Tags.DbTable));
                            assertEquals(c.getName(), c.getTag(Tags.DbColumn));
                        }
                )
        );
    }

    static Stream<SqlCase> subqueryQueries() {
        return Stream.of(

                new SqlCase(
                        "select t.col1 from (" +
                                "  select col1 from table1" +
                                ") t",
                        "Postgres",
                        new String[]{"col1"},
                        c -> {
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertNull(c.getTag(Tags.DbSchema));
                            assertEquals("col1", c.getTag(Tags.DbColumn));
                        }
                ),

                new SqlCase(
                        "select * from (" +
                                "  select col1 from table1" +
                                ") t",
                        "Postgres",
                        new String[]{"col1"},
                        c -> {
                            // Star over subquery → treated as ambiguous
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertEquals("col1", c.getTag(Tags.DbColumn));
                        }
                )
        );
    }

    static Stream<SqlCase> cteQueries() {
        return Stream.of(

                new SqlCase(
                        "with cte as (" +
                                "  select col1 from table1" +
                                ") select col1 from cte",
                        "Postgres",
                        new String[]{"col1"},
                        c -> {
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertNull(c.getTag(Tags.DbSchema));
                            assertEquals("col1", c.getTag(Tags.DbColumn));
                            assertNull(c.getTag(Tags.IsDerived));
                            assertNull(c.getTag(Tags.DbExpression));
                        }
                ),

                new SqlCase(
                        "with cte as (" +
                                "  select t.col1 from table1 t" +
                                ") select col1 from cte",
                        "Postgres",
                        new String[]{"col1"},
                        c -> {
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertNull(c.getTag(Tags.DbSchema));
                            assertEquals("col1", c.getTag(Tags.DbColumn));
                        }
                ),

                new SqlCase(
                        "with cte as (" +
                                "  select col1 as x from table1" +
                                ") select x from cte",
                        "Postgres",
                        new String[]{"x"},
                        c -> {
                            assertEquals("table1", c.getTag(Tags.DbTable));
                            assertNull(c.getTag(Tags.DbSchema));
                            assertEquals("col1", c.getTag(Tags.DbColumn));
                        }
                )
        );
    }


    @ParameterizedTest
    @MethodSource("starQueries")
    void annotatesStarQueries(SqlCase tc) {
        runCase(tc);
    }

    @ParameterizedTest
    @MethodSource("explicitQueries")
    void annotatesExplicitQueries(SqlCase tc) {
        runCase(tc);
    }


    @ParameterizedTest
    @MethodSource("derivedQueries")
    void annotatesDerivedQueries(SqlCase tc) {
        runCase(tc);
    }

    @ParameterizedTest
    @MethodSource("joinQueries")
    void annotatesJoinQueries(SqlCase tc) {
        runCase(tc);
    }

    @ParameterizedTest
    @MethodSource("quotedIdentifierQueries")
    void annotatesQuotedIdentifiers(SqlCase tc) {
        runCase(tc);
    }

    @ParameterizedTest
    @MethodSource("mssqlQueries")
    void annotatesMssqlQueries(SqlCase tc) {
        runCase(tc);
    }

    @ParameterizedTest
    @MethodSource("subqueryQueries")
    void annotatesSubqueryQueries(SqlCase tc) {
        runCase(tc);
    }

    @ParameterizedTest
    @MethodSource("cteQueries")
    void annotatesCteQueries(SqlCase tc) {
        runCase(tc);
    }


    private void runCase(SqlCase tc) {
        DataQuery q = query(tc.sql, tc.dataSource);
        DataFrame df = df(tc.columns);

        SqlAnnotator.annotate(q, df);

        for (serialization.Column<?> c : df.getColumns())
            tc.assertion.assertColumn(c);
    }
}
