package grok_connect.providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Provider;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationRunner;
import grok_connect.table_mutation.TableMutation;
import grok_connect.table_mutation.UpsertRows;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Containerized upsert tests for MySQL (connector-writes WO-4): INSERT ... ON DUPLICATE KEY UPDATE.
 * Covers the seed/merge happy path with a unique index AND the documented caveat that, without a
 * unique index, MySQL silently inserts duplicates instead of updating (no error).
 */
class MySqlMutationTest extends ContainerizedProviderBaseTest {
    protected MySqlMutationTest() {
        super(Provider.MYSQL);
    }

    private void execDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             Statement statement = c.createStatement()) {
            statement.execute(sql);
        }
    }

    private List<Object[]> queryDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             PreparedStatement statement = c.prepareStatement(sql);
             ResultSet rs = statement.executeQuery()) {
            List<Object[]> rows = new ArrayList<>();
            int columnCount = rs.getMetaData().getColumnCount();
            while (rs.next()) {
                Object[] row = new Object[columnCount];
                for (int i = 0; i < columnCount; i++)
                    row[i] = rs.getObject(i + 1);
                rows.add(row);
            }
            return rows;
        }
    }

    private long countDirect(String where) throws SQLException {
        return ((Number) queryDirect("SELECT count(*) FROM " + where).get(0)[0]).longValue();
    }

    @AfterAll
    public void dropScratchTables() throws SQLException {
        for (String table : new String[] {"mut_upsert", "mut_dup"})
            execDirect("DROP TABLE IF EXISTS " + table);
    }

    private void runDdl(String query) throws Exception {
        FuncCall call = FuncCallBuilder.fromQuery(query);
        call.func.connection = connection;
        provider.execute(call);
    }

    private MutationResult runMutation(TableMutation mutation) throws Exception {
        mutation.connection = connection;
        FuncCall call = new FuncCall();
        call.func = mutation;
        call.options = new HashMap<>();
        return MutationRunner.execute(provider, call);
    }

    @DisplayName("Upsert: seed 3, merge a 4-row payload (2 update + 2 insert) against a unique key")
    @Test
    public void upsert_seedAndMerge() throws Exception {
        runDdl("CREATE TABLE mut_upsert (id int PRIMARY KEY, region varchar(32), amount double)");
        execDirect("INSERT INTO mut_upsert VALUES (1, 'east', 10), (2, 'west', 20), (3, 'north', 30)");
        UpsertRows m = new UpsertRows();
        m.tableName = "mut_upsert";
        m.columns = Arrays.asList("id", "region", "amount");
        m.columnTypes = Arrays.asList("int", "string", "double");
        m.matchKeys = Arrays.asList("id");
        m.rows = Arrays.asList(
                Arrays.asList((Object) 1.0d, "east", 111.0d),   // update
                Arrays.asList((Object) 2.0d, "west", 222.0d),   // update
                Arrays.asList((Object) 4.0d, "south", 40.0d),   // insert
                Arrays.asList((Object) 5.0d, "central", 50.0d));// insert
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage);
        // MySQL reports 1 per inserted row + 2 per updated row (2 + 2*2 = 6); raw count is not normalized
        System.out.println("MySQL upsert affectedRows=" + result.affectedRows);
        Assertions.assertTrue(result.affectedRows >= 4, "affectedRows=" + result.affectedRows);
        Assertions.assertEquals(5, countDirect("mut_upsert"));
        Assertions.assertEquals(111.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 1").get(0)[0]).doubleValue());
        Assertions.assertEquals(222.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 2").get(0)[0]).doubleValue());
        Assertions.assertEquals(30.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 3").get(0)[0]).doubleValue()); // untouched
        Assertions.assertEquals("south", queryDirect("SELECT region FROM mut_upsert WHERE id = 4").get(0)[0]);
    }

    @DisplayName("Upsert without a unique index: MySQL silently inserts duplicates, no error (documented caveat)")
    @Test
    public void upsert_noUniqueIndex_insertsDuplicates() throws Exception {
        runDdl("CREATE TABLE mut_dup (id int, amount double)"); // NO primary/unique key
        execDirect("INSERT INTO mut_dup VALUES (1, 10)");
        UpsertRows m = new UpsertRows();
        m.tableName = "mut_dup";
        m.columns = Arrays.asList("id", "amount");
        m.columnTypes = Arrays.asList("int", "double");
        m.matchKeys = Arrays.asList("id");
        m.rows = Arrays.asList(Arrays.asList((Object) 1.0d, 99.0d));
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage); // no error even though nothing was matched
        Assertions.assertEquals(2, countDirect("mut_dup WHERE id = 1")); // duplicate inserted, not updated
    }
}
