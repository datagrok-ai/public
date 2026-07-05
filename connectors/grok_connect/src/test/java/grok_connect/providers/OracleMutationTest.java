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
 * Containerized upsert tests for Oracle (connector-writes WO-4): MERGE ... USING (SELECT ? FROM dual),
 * one row per statement executed via addBatch. Uses upper-case identifiers so the quoted mutation
 * identifiers line up with Oracle's default upper-case folding of the seed/verify SQL.
 */
class OracleMutationTest extends ContainerizedProviderBaseTest {
    protected OracleMutationTest() {
        super(Provider.ORACLE);
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
        try {
            execDirect("DROP TABLE MUT_UPSERT");
        } catch (SQLException ignore) { /* table may not exist */ }
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

    @DisplayName("Upsert: seed 3, MERGE a 4-row payload (2 update + 2 insert) by matchKeys")
    @Test
    public void upsert_seedAndMerge() throws Exception {
        try {
            execDirect("DROP TABLE MUT_UPSERT");
        } catch (SQLException ignore) { /* first run */ }
        runDdl("CREATE TABLE MUT_UPSERT (ID NUMBER PRIMARY KEY, REGION VARCHAR2(32), AMOUNT NUMBER)");
        execDirect("INSERT INTO MUT_UPSERT VALUES (1, 'east', 10)");
        execDirect("INSERT INTO MUT_UPSERT VALUES (2, 'west', 20)");
        execDirect("INSERT INTO MUT_UPSERT VALUES (3, 'north', 30)");
        UpsertRows m = new UpsertRows();
        m.tableName = "MUT_UPSERT";
        m.columns = Arrays.asList("ID", "REGION", "AMOUNT");
        m.columnTypes = Arrays.asList("int", "string", "double");
        m.matchKeys = Arrays.asList("ID");
        m.rows = Arrays.asList(
                Arrays.asList((Object) 1.0d, "east", 111.0d),   // update
                Arrays.asList((Object) 2.0d, "west", 222.0d),   // update
                Arrays.asList((Object) 4.0d, "south", 40.0d),   // insert
                Arrays.asList((Object) 5.0d, "central", 50.0d));// insert
        MutationResult result = runMutation(m);
        Assertions.assertNull(result.errorMessage);
        System.out.println("Oracle upsert affectedRows=" + result.affectedRows);
        Assertions.assertEquals(4, result.affectedRows); // 4 single-row MERGEs (SUCCESS_NO_INFO tolerated as 1)
        Assertions.assertEquals(5, countDirect("MUT_UPSERT"));
        Assertions.assertEquals(111.0d, ((Number) queryDirect("SELECT AMOUNT FROM MUT_UPSERT WHERE ID = 1").get(0)[0]).doubleValue());
        Assertions.assertEquals(222.0d, ((Number) queryDirect("SELECT AMOUNT FROM MUT_UPSERT WHERE ID = 2").get(0)[0]).doubleValue());
        Assertions.assertEquals(30.0d, ((Number) queryDirect("SELECT AMOUNT FROM MUT_UPSERT WHERE ID = 3").get(0)[0]).doubleValue()); // untouched
        Assertions.assertEquals("south", queryDirect("SELECT REGION FROM MUT_UPSERT WHERE ID = 4").get(0)[0]);
    }
}
