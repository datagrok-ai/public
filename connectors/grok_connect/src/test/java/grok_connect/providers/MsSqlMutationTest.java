package grok_connect.providers;

import grok_connect.connectors_info.Credentials;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Provider;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationRunner;
import grok_connect.table_mutation.TableMutation;
import grok_connect.table_mutation.UpsertRows;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
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
 * Containerized upsert tests for MS SQL (connector-writes WO-4): chunked MERGE-over-VALUES with the
 * 2100-parameter limit handled by {@link grok_connect.providers.MsSqlDataProvider#upsertBatchRows}.
 */
class MsSqlMutationTest extends ContainerizedProviderBaseTest {
    private static final String DEFAULT_DATABASE_NAME = "master"; // can't read the db name off the MSSQL container

    protected MsSqlMutationTest() {
        super(Provider.MSSQL);
    }

    @Override
    @BeforeEach
    public void beforeEach() {
        credentials = new Credentials();
        credentials.parameters.put(DbCredentials.LOGIN, container.getUsername());
        credentials.parameters.put(DbCredentials.PASSWORD, container.getPassword());
        connection = new DataConnection();
        connection.credentials = credentials;
        connection.dataSource = provider.descriptor.type;
        connection.parameters.put(DbCredentials.SERVER, container.getHost());
        connection.parameters.put(DbCredentials.PORT, (double) container.getFirstMappedPort());
        connection.parameters.put(DbCredentials.DB, DEFAULT_DATABASE_NAME);
    }

    private String url() {
        return "jdbc:sqlserver://" + container.getHost() + ":" + container.getFirstMappedPort()
                + ";databaseName=" + DEFAULT_DATABASE_NAME + ";encrypt=false;trustServerCertificate=true";
    }

    private void execDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(url(), container.getUsername(), container.getPassword());
             Statement statement = c.createStatement()) {
            statement.execute(sql);
        }
    }

    private List<Object[]> queryDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(url(), container.getUsername(), container.getPassword());
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
        execDirect("IF OBJECT_ID('mut_upsert', 'U') IS NOT NULL DROP TABLE mut_upsert");
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
        execDirect("IF OBJECT_ID('mut_upsert', 'U') IS NOT NULL DROP TABLE mut_upsert");
        runDdl("CREATE TABLE mut_upsert (id int PRIMARY KEY, region nvarchar(32), amount float)");
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
        System.out.println("MS SQL upsert affectedRows=" + result.affectedRows);
        Assertions.assertEquals(4, result.affectedRows); // MERGE reports total rows touched
        Assertions.assertEquals(5, countDirect("mut_upsert"));
        Assertions.assertEquals(111.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 1").get(0)[0]).doubleValue());
        Assertions.assertEquals(222.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 2").get(0)[0]).doubleValue());
        Assertions.assertEquals(30.0d, ((Number) queryDirect("SELECT amount FROM mut_upsert WHERE id = 3").get(0)[0]).doubleValue()); // untouched
        Assertions.assertEquals("south", queryDirect("SELECT region FROM mut_upsert WHERE id = 4").get(0)[0]);
    }
}
