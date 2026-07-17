package grok_connect.providers;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Provider;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryManager;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import serialization.DataFrame;

/**
 * Integration tests for the §6.2 raw-SQL policy (connector-writes WO-B13) against a containerized
 * Postgres. Layered enforcement is pinned explicitly: the first-keyword classifier refuses obvious
 * writes with the structured {@code {error: 'read-only'}} shape, and the driver read-only session
 * stops what the classifier cannot see (a write hidden in a CTE). With the {@code auditRawWrites}
 * option (allowRawWrites ON — today's behavior), writes still execute but are detected post-hoc.
 */
class ReadOnlyQueryTest extends ContainerizedProviderBaseTest {

    protected ReadOnlyQueryTest() {
        super(Provider.POSTGRESQL);
    }

    @BeforeAll
    public void initEnvironment() throws SQLException {
        // QueryManager resolves its provider through the GrokConnect static (streaming-path test)
        GrokConnect.providerManager = new ProviderManager();
        execDirect("CREATE TABLE IF NOT EXISTS ro_probe (id int)");
    }

    @AfterAll
    public void dropScratchTables() throws SQLException {
        execDirect("DROP TABLE IF EXISTS ro_probe");
    }

    /** Direct JDBC, bypassing grok_connect — verification must not depend on the code under test. */
    private void execDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             Statement statement = c.createStatement()) {
            statement.execute(sql);
        }
    }

    private long countProbe() throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             ResultSet rs = c.createStatement().executeQuery("SELECT count(*) FROM ro_probe")) {
            rs.next();
            return rs.getLong(1);
        }
    }

    private FuncCall queryCall(String sql, String option) {
        FuncCall call = FuncCallBuilder.fromQuery(sql);
        call.func.connection = connection;
        if (option != null)
            call.options.put(option, Boolean.TRUE);
        return call;
    }

    private String refusalMessage(GrokConnectException ex) {
        return ex.getCause() != null ? ex.getCause().getMessage() : ex.getMessage();
    }

    @DisplayName("readOnly: SELECT executes normally")
    @Test
    public void readOnly_selectWorks() throws Exception {
        DataFrame result = provider.execute(queryCall("SELECT 1 AS one", DataProvider.READ_ONLY));
        Assertions.assertEquals(1, result.rowCount);
    }

    @DisplayName("readOnly: a SELECT with a leading comment classifies as a read")
    @Test
    public void readOnly_selectWithLeadingComment() throws Exception {
        DataFrame result = provider.execute(queryCall("-- reporting query\nSELECT 1 AS one", DataProvider.READ_ONLY));
        Assertions.assertEquals(1, result.rowCount);
    }

    @DisplayName("readOnly: INSERT is refused by the classifier with the structured read-only shape")
    @Test
    public void readOnly_insertRefusedByClassifier() throws Exception {
        long before = countProbe();
        GrokConnectException ex = Assertions.assertThrows(GrokConnectException.class,
                () -> provider.execute(queryCall("INSERT INTO ro_probe (id) VALUES (1)", DataProvider.READ_ONLY)));
        String message = refusalMessage(ex);
        Assertions.assertTrue(message.contains("\"error\":\"read-only\""),
                "expected the structured refusal shape, got: " + message);
        Assertions.assertTrue(message.contains("allowRawWrites"), "the refusal must name the flag: " + message);
        Assertions.assertTrue(message.contains("DataConnection.AddRows"),
                "the refusal must name the fine privileges: " + message);
        Assertions.assertEquals(before, countProbe(), "no row must be written");
    }

    @DisplayName("readOnly: batch mode classifies every statement, not just the first")
    @Test
    public void readOnly_batchModeChecksEveryStatement() throws Exception {
        long before = countProbe();
        FuncCall call = queryCall("SELECT 1 AS one\n--batch\nINSERT INTO ro_probe (id) VALUES (2)",
                DataProvider.READ_ONLY);
        call.func.options = new HashMap<>();
        call.func.options.put("batchMode", "true");
        GrokConnectException ex = Assertions.assertThrows(GrokConnectException.class, () -> provider.execute(call));
        Assertions.assertTrue(refusalMessage(ex).contains("\"error\":\"read-only\""),
                "the second (write) statement of the batch must be refused, got: " + refusalMessage(ex));
        Assertions.assertEquals(before, countProbe(), "no row must be written");
    }

    @DisplayName("readOnly: a CTE-hidden INSERT passes the classifier but is stopped by the read-only session (layered proof)")
    @Test
    public void readOnly_cteBypassStoppedBySession() throws Exception {
        long before = countProbe();
        GrokConnectException ex = Assertions.assertThrows(GrokConnectException.class,
                () -> provider.execute(queryCall(
                        "WITH ins AS (INSERT INTO ro_probe (id) VALUES (42) RETURNING id) SELECT * FROM ins",
                        DataProvider.READ_ONLY)));
        String message = refusalMessage(ex);
        // NOT the classifier's structured shape — the refusal comes from the database session (layer 2)
        Assertions.assertFalse(message.contains("\"error\":\"read-only\""),
                "the CTE must pass the classifier; the session must refuse it, got: " + message);
        Assertions.assertTrue(message.toLowerCase().contains("read-only"),
                "expected the PG read-only-transaction error, got: " + message);
        Assertions.assertEquals(before, countProbe(), "the CTE-hidden insert must not persist");
    }

    @DisplayName("readOnly: the streaming path (QueryManager) refuses an INSERT at initResultSet")
    @Test
    public void readOnly_streamingPathRefused() throws Exception {
        long before = countProbe();
        FuncCall call = queryCall("INSERT INTO ro_probe (id) VALUES (3)", DataProvider.READ_ONLY);
        call.func.type = "DataQuery";
        QueryManager queryManager = new QueryManager(GrokConnect.gson.toJson(call));
        SQLException ex = Assertions.assertThrows(SQLException.class,
                () -> queryManager.initResultSet(queryManager.getQuery()));
        Assertions.assertTrue(ex.getMessage().contains("\"error\":\"read-only\""),
                "expected the structured refusal on the streaming path, got: " + ex.getMessage());
        queryManager.close(false);
        Assertions.assertEquals(before, countProbe(), "no row must be written");
    }

    @DisplayName("auditRawWrites (flag ON): INSERT executes AND is detected post-hoc")
    @Test
    public void auditFlag_insertExecutesAndIsDetected() throws Exception {
        execDirect("DELETE FROM ro_probe");
        FuncCall call = queryCall("INSERT INTO ro_probe (id) VALUES (7)", DataProvider.AUDIT_RAW_WRITES);
        provider.execute(call);
        Assertions.assertEquals(1, countProbe(), "flag ON preserves today's behavior - the write lands");
        Assertions.assertEquals(Boolean.TRUE, call.aux.get(DataProvider.RAW_WRITE_DETECTED),
                "the raw write must be detected for the audit event");
        execDirect("DELETE FROM ro_probe");
    }

    @DisplayName("auditRawWrites (flag ON): DDL is detected too")
    @Test
    public void auditFlag_ddlDetected() throws Exception {
        FuncCall call = queryCall("CREATE TABLE ro_audit_ddl (id int)", DataProvider.AUDIT_RAW_WRITES);
        try {
            provider.execute(call);
            Assertions.assertEquals(Boolean.TRUE, call.aux.get(DataProvider.RAW_WRITE_DETECTED));
        } finally {
            execDirect("DROP TABLE IF EXISTS ro_audit_ddl");
        }
    }

    @DisplayName("auditRawWrites (flag ON): the streaming path detects and exposes the flag")
    @Test
    public void auditFlag_streamingPathDetects() throws Exception {
        execDirect("DELETE FROM ro_probe");
        FuncCall call = queryCall("INSERT INTO ro_probe (id) VALUES (8)", DataProvider.AUDIT_RAW_WRITES);
        call.func.type = "DataQuery";
        QueryManager queryManager = new QueryManager(GrokConnect.gson.toJson(call));
        queryManager.initResultSet(queryManager.getQuery());
        Assertions.assertTrue(queryManager.isRawWriteDetected(),
                "QueryManager must expose the detection for the COMPLETED RAW_WRITE token");
        queryManager.close(true);
        Assertions.assertEquals(1, countProbe(), "the streamed insert commits on close(true)");
        execDirect("DELETE FROM ro_probe");
    }

    @DisplayName("auditRawWrites (flag ON): a SELECT is not flagged")
    @Test
    public void auditFlag_selectNotDetected() throws Exception {
        FuncCall call = queryCall("SELECT 1 AS one", DataProvider.AUDIT_RAW_WRITES);
        provider.execute(call);
        Assertions.assertNull(call.aux.get(DataProvider.RAW_WRITE_DETECTED));
    }

    @DisplayName("no options: raw INSERT works and nothing is flagged (today's behavior preserved)")
    @Test
    public void noOptions_insertWorksUndetected() throws Exception {
        execDirect("DELETE FROM ro_probe");
        FuncCall call = queryCall("INSERT INTO ro_probe (id) VALUES (9)", null);
        provider.execute(call);
        Assertions.assertEquals(1, countProbe());
        Assertions.assertNull(call.aux.get(DataProvider.RAW_WRITE_DETECTED));
        execDirect("DELETE FROM ro_probe");
    }

    @DisplayName("descriptor honesty: Postgres declares real read-only session enforcement")
    @Test
    public void descriptor_declaresEnforcement() {
        Assertions.assertTrue(provider.descriptor.readOnlySessionEnforced);
    }
}
