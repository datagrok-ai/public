package grok_connect.handlers;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Arrays;
import java.util.HashMap;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.ContainerizedProviderBaseTest;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Provider;
import grok_connect.table_mutation.InsertRows;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationRunner;
import grok_connect.utils.ProviderManager;
import org.eclipse.jetty.websocket.api.RemoteEndpoint;
import org.eclipse.jetty.websocket.api.Session;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;

/**
 * GROK-20322 regression: a client that dies mid-MUTATION-stream (socket close or error without
 * {@code MUTATION EOF}) must roll the transaction back BEFORE the connection returns to the pool
 * (the phase-B e2e Hikari trap: a poisoned borrower saw stale rows and no-effect DDL). Drives the
 * REAL {@link SessionHandler} WS state machine with a mocked Jetty session against a containerized
 * Postgres, then proves fresh borrowers from the SAME pool see a clean transaction: the partial
 * bulk data is gone, an inline insert lands with correct affectedRows, and DDL takes effect.
 */
class SessionHandlerMutationAbortTest extends ContainerizedProviderBaseTest {

    protected SessionHandlerMutationAbortTest() {
        super(Provider.POSTGRESQL);
    }

    @BeforeAll
    public void initProviderManager() {
        GrokConnect.providerManager = new ProviderManager();
    }

    @AfterAll
    public void dropScratchTables() throws SQLException {
        for (String table : new String[] {"sh_abort_close", "sh_abort_error", "sh_probe_ddl"})
            execDirect("DROP TABLE IF EXISTS " + table);
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

    private boolean tableExistsDirect(String table) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             ResultSet rs = c.createStatement().executeQuery(
                     "SELECT count(*) FROM information_schema.tables WHERE table_name = '" + table + "'")) {
            rs.next();
            return rs.getLong(1) > 0;
        }
    }

    private SessionHandler newHandler() {
        Session session = Mockito.mock(Session.class);
        RemoteEndpoint remote = Mockito.mock(RemoteEndpoint.class);
        Mockito.when(session.getRemote()).thenReturn(remote);
        return new SessionHandler(session, true);
    }

    private String mutationHeader(String table) {
        InsertRows m = new InsertRows();
        m.tableName = table;
        m.columns = Arrays.asList("id", "note");
        m.columnTypes = Arrays.asList("int", "string");
        m.bulk = true;
        m.rows = null;
        m.connection = connection;
        FuncCall call = new FuncCall();
        call.func = m;
        call.options = new HashMap<>();
        return "MUTATION " + GrokConnect.gson.toJson(call);
    }

    private byte[] chunk(int from, int count) {
        Integer[] id = new Integer[count];
        String[] note = new String[count];
        for (int i = 0; i < count; i++) {
            id[i] = from + i;
            note[i] = "n" + (from + i);
        }
        return DataFrame.fromColumns(new IntColumn("id", id), new StringColumn("note", note)).toByteArray();
    }

    /** The fresh-borrower probes: the exact operations the poisoned pool silently broke live. */
    private void assertPoolIsClean(String abortedTable) throws Exception {
        // (1) partial bulk data must be gone — the abort rolled back, nothing was committed by the pool reset
        Assertions.assertEquals(0, countDirect(abortedTable), "mid-stream abort must leave zero rows");

        // (2) an inline insert on a fresh borrower reports real affectedRows and the rows are visible in psql
        InsertRows ins = new InsertRows();
        ins.tableName = abortedTable;
        ins.columns = Arrays.asList("id", "note");
        ins.columnTypes = Arrays.asList("int", "string");
        ins.rows = Arrays.asList(
                Arrays.asList((Object) 1d, "alpha"),
                Arrays.asList((Object) 2d, "beta"));
        ins.connection = connection;
        FuncCall call = new FuncCall();
        call.func = ins;
        call.options = new HashMap<>();
        MutationResult result = MutationRunner.execute(provider, call);
        Assertions.assertNull(result.errorMessage, "fresh-borrower insert must not fail");
        Assertions.assertEquals(2, result.affectedRows, "fresh-borrower insert must report real affectedRows");
        Assertions.assertEquals(2, countDirect(abortedTable), "fresh-borrower insert must be committed and visible");

        // (3) DDL on fresh borrowers takes effect (create, then drop, each verified via direct JDBC)
        FuncCall create = FuncCallBuilder.fromQuery("CREATE TABLE sh_probe_ddl (id int)");
        create.func.connection = connection;
        provider.execute(create);
        Assertions.assertTrue(tableExistsDirect("sh_probe_ddl"), "DDL create on a fresh borrower must take effect");
        FuncCall drop = FuncCallBuilder.fromQuery("DROP TABLE sh_probe_ddl");
        drop.func.connection = connection;
        provider.execute(drop);
        Assertions.assertFalse(tableExistsDirect("sh_probe_ddl"), "DDL drop on a fresh borrower must take effect");
    }

    @DisplayName("Mid-stream socket close (no EOF): the session teardown rolls back before pool return; fresh borrowers are clean")
    @Test
    public void socketCloseMidStream_rollsBackAndPoolStaysClean() throws Throwable {
        execDirect("DROP TABLE IF EXISTS sh_abort_close");
        execDirect("CREATE TABLE sh_abort_close (id int PRIMARY KEY, note text)");
        SessionHandler handler = newHandler();
        handler.onMessage(mutationHeader("sh_abort_close"));
        byte[] part = chunk(0, 500);
        handler.onBinary(part, 0, part.length); // rows now sit in the open COPY transaction
        handler.onClose(); // the client died — Jetty delivers the close without MUTATION EOF
        assertPoolIsClean("sh_abort_close");
    }

    @DisplayName("Mid-stream socket error: onError tears the mutation down with rollback; fresh borrowers are clean")
    @Test
    public void socketErrorMidStream_rollsBackAndPoolStaysClean() throws Throwable {
        execDirect("DROP TABLE IF EXISTS sh_abort_error");
        execDirect("CREATE TABLE sh_abort_error (id int PRIMARY KEY, note text)");
        SessionHandler handler = newHandler();
        handler.onMessage(mutationHeader("sh_abort_error"));
        byte[] part = chunk(0, 500);
        handler.onBinary(part, 0, part.length);
        handler.onError(new java.util.concurrent.TimeoutException("Idle timeout expired")); // the idle-timeout path
        handler.onClose(); // Jetty always follows onError with onClose
        assertPoolIsClean("sh_abort_error");
    }
}
