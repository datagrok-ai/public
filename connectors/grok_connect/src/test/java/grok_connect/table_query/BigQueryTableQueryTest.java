package grok_connect.table_query;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.providers.BigQueryDataProvider;
import grok_connect.utils.ProviderManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;

import java.util.ArrayList;
import java.util.Arrays;

@TestInstance(TestInstance.Lifecycle.PER_CLASS)
public class BigQueryTableQueryTest {
    private static final String PROJECT = "my-gall-1996";
    private static final String DATASET = "datagrok_test";
    private static final String OTHER_DATASET = "other_ds";
    private static final String TABLE = "people";

    private final BigQueryDataProvider provider = new BigQueryDataProvider();
    private TableQuery query;

    @BeforeAll
    public void setUpAll() {
        GrokConnect.providerManager = new ProviderManager();
    }

    @BeforeEach
    public void setUp() {
        query = new TableQuery();
        query.tableName = TABLE;
        query.schema = DATASET;
        query.params = new ArrayList<>();
        query.connection = new DataConnection();
        query.connection.dataSource = provider.descriptor.type;
        query.connection.parameters.put("projectId", PROJECT);
    }

    @Test
    @DisplayName("FROM is rendered as a single backticked `project.dataset.table`")
    public void mainTableIsProjectQualified() {
        String sql = provider.queryTableSql(query);
        String fullPath = "`" + PROJECT + "." + DATASET + "." + TABLE + "`";
        Assertions.assertTrue(
                sql.contains("FROM" + System.lineSeparator() + fullPath),
                "FROM should be a single backticked project-qualified path:\n" + sql);
        Assertions.assertFalse(sql.contains("`" + PROJECT + "`.`"),
                "segments should not be split-bracketed (BQ rejects `project-id` with dashes):\n" + sql);
    }

    @Test
    @DisplayName("Fields shaped like `dataset.table.col` gain the project prefix and merge into one identifier")
    public void threePartFieldsGetProjectPrefix() {
        query.fields = Arrays.asList(
                DATASET + "." + TABLE + ".id",
                DATASET + "." + TABLE + ".name");
        String sql = provider.queryTableSql(query);
        String tablePath = "`" + PROJECT + "." + DATASET + "." + TABLE + "`";
        Assertions.assertTrue(sql.contains(tablePath + ".`id`"),
                "column id should reference the merged table path:\n" + sql);
        Assertions.assertTrue(sql.contains(tablePath + ".`name`"),
                "column name should reference the merged table path:\n" + sql);
    }

    @Test
    @DisplayName("Cross-dataset JOIN picks up the joined dataset and prefixes it too")
    public void crossDatasetJoinGetsProjectPrefix() {
        TableJoin join = new TableJoin();
        join.joinType = "left";
        join.leftTableName = TABLE;
        join.rightTableName = OTHER_DATASET + ".orders";
        join.leftTableKeys = new ArrayList<String>() {{ add("id"); }};
        join.rightTableKeys = new ArrayList<String>() {{ add("person_id"); }};
        query.joins.add(join);

        String sql = provider.queryTableSql(query);
        Assertions.assertTrue(sql.contains("`" + PROJECT + "." + OTHER_DATASET + ".orders`"),
                "joined table in other dataset should be a single backticked path:\n" + sql);
        Assertions.assertTrue(sql.contains("`" + PROJECT + "." + DATASET + "." + TABLE + "`"),
                "main table should still be a single backticked path:\n" + sql);
    }

    @Test
    @DisplayName("Paths already starting with the project are not double-prefixed")
    public void alreadyQualifiedPathsAreNotDoubled() {
        query.schema = PROJECT + "." + DATASET;
        String sql = provider.queryTableSql(query);
        Assertions.assertFalse(sql.contains(PROJECT + "." + PROJECT),
                "project prefix should not appear twice:\n" + sql);
        Assertions.assertTrue(sql.contains("`" + PROJECT + "." + DATASET + "." + TABLE + "`"),
                "FROM should be a single backticked project-qualified path:\n" + sql);
    }

    @Test
    @DisplayName("No projectId on connection leaves the SQL untouched")
    public void missingProjectIdIsNoOp() {
        query.connection.parameters.remove("projectId");
        String sql = provider.queryTableSql(query);
        Assertions.assertFalse(sql.contains("`" + PROJECT + "`") || sql.contains(PROJECT),
                "no project prefix should be injected when projectId is absent:\n" + sql);
        Assertions.assertTrue(sql.contains("`" + DATASET + "`.`" + TABLE + "`"),
                "split-bracketed dataset-qualified FROM should remain (non-BQ default bracketing):\n" + sql);
    }
}
