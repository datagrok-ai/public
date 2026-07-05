package grok_connect.providers;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Provider;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import serialization.DataFrame;

import java.sql.SQLException;
import java.util.HashMap;

/**
 * Regression tests for GROK-20323: one execution = one transaction with explicit
 * commit-on-success / rollback-on-error on both query paths (non-streaming
 * JdbcDataProvider.execute and streaming QueryManager).
 */
class TransactionSemanticsTest extends ContainerizedProviderBaseTest {
    protected TransactionSemanticsTest() {
        super(Provider.POSTGRESQL);
    }

    @BeforeAll
    public void initGrokConnectStatics() {
        // QueryManager resolves its provider through GrokConnect.providerManager,
        // which is only set in GrokConnect.main
        if (GrokConnect.providerManager == null)
            GrokConnect.providerManager = new ProviderManager();
    }

    private DataFrame run(String query) throws Exception {
        FuncCall call = FuncCallBuilder.fromQuery(query);
        call.func.connection = connection;
        return provider.execute(call);
    }

    private FuncCall buildBatchCall(String... queries) {
        FuncCall call = FuncCallBuilder.fromQuery(String.join("\n--batch\n", queries));
        call.func.connection = connection;
        call.func.options = new HashMap<>();
        call.func.options.put("batchMode", "true");
        return call;
    }

    private QueryManager buildQueryManager(FuncCall call) {
        call.func.type = "DataQuery"; // DataQueryDeserializer dispatches on #type
        return new QueryManager(GrokConnect.gson.toJson(call));
    }

    @DisplayName("(a) Non-streaming execute commits: DDL + INSERT persist across connections")
    @Test
    public void execute_write_isCommitted() throws Exception {
        run("CREATE TABLE tx_exec_ok (id int PRIMARY KEY)");
        run("INSERT INTO tx_exec_ok (id) VALUES (1)");
        DataFrame df = run("SELECT id FROM tx_exec_ok");
        Assertions.assertEquals(1, df.rowCount);
    }

    @DisplayName("(b) Streaming QueryManager.close(true) commits an INSERT")
    @Test
    public void queryManager_closeWithCommit_persistsWrite() throws Exception {
        run("CREATE TABLE tx_qm_ok (id int PRIMARY KEY)");
        FuncCall call = FuncCallBuilder.fromQuery("INSERT INTO tx_qm_ok (id) VALUES (1)");
        call.func.connection = connection;
        QueryManager qm = buildQueryManager(call);
        try {
            qm.initResultSet(qm.getQuery());
            if (qm.isResultSetInitialized())
                qm.getSubDF(1);
            qm.close(true);
        } catch (Throwable t) {
            qm.close(false);
            throw t;
        }
        DataFrame df = run("SELECT id FROM tx_qm_ok");
        Assertions.assertEquals(1, df.rowCount);
    }

    @DisplayName("(c1) Streaming error path: close(false) after a mid-batch failure rolls back the first statement")
    @Test
    public void queryManager_closeAfterMidBatchFailure_rollsBack() throws Exception {
        run("CREATE TABLE tx_qm_err (id int PRIMARY KEY)");
        FuncCall call = buildBatchCall(
                "INSERT INTO tx_qm_err (id) VALUES (1)",
                "INSERT INTO tx_qm_err (id) VALUES (1)"); // PK violation
        QueryManager qm = buildQueryManager(call);
        Assertions.assertThrows(SQLException.class, () -> qm.initResultSet(qm.getQuery()));
        qm.close(false);
        DataFrame df = run("SELECT id FROM tx_qm_err");
        Assertions.assertEquals(0, df.rowCount);
    }

    @DisplayName("(c2) Streaming abort before COMPLETED_OK: close(false) on a healthy transaction rolls back, not commits")
    @Test
    public void queryManager_closeWithoutCompletion_rollsBack() throws Exception {
        run("CREATE TABLE tx_qm_abort (id int PRIMARY KEY)");
        FuncCall call = buildBatchCall(
                "INSERT INTO tx_qm_abort (id) VALUES (1)",
                "SELECT id FROM tx_qm_abort");
        QueryManager qm = buildQueryManager(call);
        qm.initResultSet(qm.getQuery());
        // simulate SessionHandler error/abort path: session torn down before COMPLETED_OK
        qm.close(false);
        DataFrame df = run("SELECT id FROM tx_qm_abort");
        Assertions.assertEquals(0, df.rowCount);
    }

    @DisplayName("(d) Non-streaming execute rolls back on a mid-batch failure")
    @Test
    public void execute_midBatchFailure_rollsBack() throws Exception {
        run("CREATE TABLE tx_exec_err (id int PRIMARY KEY)");
        FuncCall call = buildBatchCall(
                "INSERT INTO tx_exec_err (id) VALUES (1)",
                "INSERT INTO tx_exec_err (id) VALUES (1)"); // PK violation
        Assertions.assertThrows(GrokConnectException.class, () -> provider.execute(call));
        DataFrame df = run("SELECT id FROM tx_exec_err");
        Assertions.assertEquals(0, df.rowCount);
    }
}
