package grok_connect.handlers;

import java.lang.reflect.Field;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.ContainerizedProviderBaseTest;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Provider;
import grok_connect.utils.ProviderManager;
import org.eclipse.jetty.websocket.api.RemoteEndpoint;
import org.eclipse.jetty.websocket.api.Session;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;

/**
 * GROK-20322 (WO-B14) regression: a streaming adhoc query's effects must be committed BEFORE the
 * {@code COMPLETED} token is sent — Datlas resolves the caller's future on that token, so the old
 * commit-at-WS-close contract lost a ~4 ms visibility race to follow-up operations on other pooled
 * connections. Drives the REAL {@link SessionHandler} WS state machine with a mocked Jetty session
 * against a containerized Postgres; the mocked remote runs a fresh-connection probe at the exact
 * moment the {@code COMPLETED} send is enqueued, which is the happens-before edge under test.
 */
class SessionHandlerQueryCommitVisibilityTest extends ContainerizedProviderBaseTest {
    private List<String> sent;
    private Callable<Long> completedProbe;
    private Long observedAtCompleted;
    private SessionHandler lastHandler;

    protected SessionHandlerQueryCommitVisibilityTest() {
        super(Provider.POSTGRESQL);
    }

    @BeforeAll
    public void initGrokConnectStatics() throws Exception {
        GrokConnect.providerManager = new ProviderManager();
        // the SIZE RECEIVED branch submits the next-chunk future through GrokConnect's pool,
        // which only GrokConnect.main creates
        Field poolField = GrokConnect.class.getDeclaredField("threadPool");
        poolField.setAccessible(true);
        if (poolField.get(null) == null)
            poolField.set(null, Executors.newCachedThreadPool());
    }

    @AfterAll
    public void dropScratchTables() throws SQLException {
        for (String table : new String[] {"wo14_ins", "wo14_ddl", "wo14_rollback", "wo14_defer"})
            execDirect("DROP TABLE IF EXISTS " + table);
    }

    @BeforeEach
    public void resetObservation() {
        sent = new ArrayList<>();
        completedProbe = null;
        observedAtCompleted = null;
    }

    @AfterEach
    public void closeHandlerQuietly() {
        // a failed assertion leaves the handler's pooled connection mid-transaction; close it so
        // its locks cannot block dropScratchTables (an unfixed-code run would hang on DROP TABLE)
        if (lastHandler != null) {
            try {
                lastHandler.onClose();
            } catch (Throwable ignore) {}
            lastHandler = null;
        }
    }

    private void execDirect(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             Statement statement = c.createStatement()) {
            statement.execute(sql);
        }
    }

    private long countDirect(String table) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             ResultSet rs = c.createStatement().executeQuery("SELECT count(*) FROM " + table)) {
            rs.next();
            return rs.getLong(1);
        }
    }

    private long tableCountDirect(String table) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             ResultSet rs = c.createStatement().executeQuery(
                     "SELECT count(*) FROM information_schema.tables WHERE table_name = '" + table + "'")) {
            rs.next();
            return rs.getLong(1);
        }
    }

    /** The mocked remote records every outbound frame and runs the visibility probe on a fresh
     *  connection the moment the COMPLETED send is enqueued — BEFORE COMPLETED_OK/EOF/onClose. */
    private SessionHandler newHandler() {
        Session session = Mockito.mock(Session.class);
        RemoteEndpoint remote = Mockito.mock(RemoteEndpoint.class);
        Mockito.when(session.getRemote()).thenReturn(remote);
        Mockito.when(remote.sendStringByFuture(Mockito.anyString())).thenAnswer(invocation -> {
            String message = invocation.getArgument(0);
            sent.add(message);
            if (message.startsWith("COMPLETED ") && completedProbe != null)
                observedAtCompleted = completedProbe.call();
            return null;
        });
        lastHandler = new SessionHandler(session, true);
        return lastHandler;
    }

    private void sendQueryHeader(SessionHandler handler, FuncCall call) throws Throwable {
        call.func.connection = connection;
        call.func.type = "DataQuery"; // DataQueryDeserializer dispatches on #type
        handler.onMessage("QUERY " + GrokConnect.gson.toJson(call));
    }

    /** Full client half of the protocol for a no-resultset statement: header, size ack, chunk ack. */
    private void driveQuery(SessionHandler handler, String sql) throws Throwable {
        sendQueryHeader(handler, FuncCallBuilder.fromQuery(sql));
        handler.onMessage("DATAFRAME PART SIZE RECEIVED");
        handler.onMessage("PART OK");
    }

    @DisplayName("INSERT adhoc: the row is committed and visible from a fresh connection the moment COMPLETED is sent; close(true) after it is an idempotent no-op")
    @Test
    public void insert_visibleAtCompletedToken_closeIdempotent() throws Throwable {
        execDirect("DROP TABLE IF EXISTS wo14_ins");
        execDirect("CREATE TABLE wo14_ins (id int PRIMARY KEY)");
        completedProbe = () -> countDirect("wo14_ins");
        SessionHandler handler = newHandler();
        driveQuery(handler, "INSERT INTO wo14_ins (id) VALUES (1)");
        Assertions.assertEquals(Long.valueOf(1), observedAtCompleted,
                "the row must be visible to a fresh connection when COMPLETED is sent (before onClose)");
        handler.onMessage("COMPLETED_OK");
        handler.onClose(); // close(true) — commit already done, must neither throw nor double-apply
        Assertions.assertEquals(1, countDirect("wo14_ins"), "the row must be present exactly once after close");
    }

    @DisplayName("CREATE TABLE adhoc (the matrix scenario): the table is visible from a fresh connection the moment COMPLETED is sent")
    @Test
    public void createTable_visibleAtCompletedToken() throws Throwable {
        execDirect("DROP TABLE IF EXISTS wo14_ddl");
        completedProbe = () -> tableCountDirect("wo14_ddl");
        SessionHandler handler = newHandler();
        driveQuery(handler, "CREATE TABLE wo14_ddl (id int)");
        Assertions.assertEquals(Long.valueOf(1), observedAtCompleted,
                "the table must be visible to a fresh connection when COMPLETED is sent (before onClose)");
        handler.onMessage("COMPLETED_OK");
        handler.onClose();
        Assertions.assertEquals(1, tableCountDirect("wo14_ddl"), "the table must still exist after close");
    }

    @DisplayName("Failing statement: ERROR sent, no COMPLETED, nothing committed (rollback-on-error preserved)")
    @Test
    public void failingStatement_rollsBackAndSendsError() throws Throwable {
        execDirect("DROP TABLE IF EXISTS wo14_rollback");
        execDirect("CREATE TABLE wo14_rollback (id int PRIMARY KEY)");
        SessionHandler handler = newHandler();
        FuncCall call = FuncCallBuilder.fromQuery(
                "INSERT INTO wo14_rollback (id) VALUES (1)\n--batch\nINSERT INTO wo14_rollback (id) VALUES (1)"); // PK violation mid-batch
        call.func.options = new HashMap<>();
        call.func.options.put("batchMode", "true");
        Throwable error = Assertions.assertThrows(SQLException.class, () -> sendQueryHeader(handler, call));
        handler.onError(error); // Jetty follows the onMessage throw with onError
        Assertions.assertTrue(sent.stream().anyMatch(m -> m.startsWith("ERROR:")), "an ERROR frame must reach the server");
        Assertions.assertTrue(sent.stream().noneMatch(m -> m.startsWith("COMPLETED ")), "no COMPLETED token after a failure");
        Assertions.assertEquals(0, countDirect("wo14_rollback"), "the successful first statement must be rolled back");
        handler.onClose(); // Jetty always follows onError with onClose — must not throw
    }

    @DisplayName("Failed commit (deferred constraint) surfaces as a query ERROR before any COMPLETED token, not swallowed post-success")
    @Test
    public void failedCommit_surfacesAsError() throws Throwable {
        execDirect("DROP TABLE IF EXISTS wo14_defer");
        execDirect("CREATE TABLE wo14_defer (id int, CONSTRAINT wo14_defer_uniq UNIQUE (id) DEFERRABLE INITIALLY DEFERRED)");
        execDirect("INSERT INTO wo14_defer (id) VALUES (1)");
        SessionHandler handler = newHandler();
        // the duplicate passes the deferred check at statement time; the violation fires at COMMIT
        sendQueryHeader(handler, FuncCallBuilder.fromQuery("INSERT INTO wo14_defer (id) VALUES (1)"));
        handler.onMessage("DATAFRAME PART SIZE RECEIVED");
        Throwable error = Assertions.assertThrows(SQLException.class, () -> handler.onMessage("PART OK"));
        handler.onError(error);
        Assertions.assertTrue(sent.stream().anyMatch(m -> m.startsWith("ERROR:")), "a failed commit must reach the server as ERROR");
        Assertions.assertTrue(sent.stream().noneMatch(m -> m.startsWith("COMPLETED ")), "no COMPLETED token when the commit failed");
        Assertions.assertEquals(1, countDirect("wo14_defer"), "only the seed row must remain");
        handler.onClose();
    }
}
