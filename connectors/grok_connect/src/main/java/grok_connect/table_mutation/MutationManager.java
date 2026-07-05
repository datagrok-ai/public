package grok_connect.table_mutation;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.HashMap;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.GrokConnectUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Drives one streamed bulk mutation over the WebSocket transport (connector-writes WO-5) — the
 * QueryManager sibling for writes. Owns its own connection, transaction and {@link BulkLoader}
 * lifecycle: {@link #start()} opens the connection and loader, {@link #feed} consumes CSV chunks,
 * {@link #finish()} flushes and commits, {@link #abort()} rolls back. The connection is always
 * closed before it returns to the pool (GROK-20323 contract). Unlike QueryManager this never
 * commits on close — commit happens only in {@link #finish()}.
 */
public class MutationManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(MutationManager.class);

    public final JdbcDataProvider provider;
    private final InsertRows mutation;
    private Connection connection;
    private BulkLoader loader;
    private boolean finished;

    public MutationManager(String message) {
        FuncCall call = GrokConnect.gson.fromJson(message, FuncCall.class);
        if (call == null || call.func == null)
            throw new MutationValidationException("Bulk mutation header must be a FuncCall with an InsertRows func");
        if (call.options == null)
            call.options = new HashMap<>();
        call.setParamValues();
        call.afterDeserialization();
        if (!(call.func instanceof InsertRows))
            throw new MutationValidationException("Bulk mutation requires an InsertRows func, got: " + call.func.type);
        mutation = (InsertRows) call.func;
        if (!mutation.bulk)
            throw new MutationValidationException("Bulk mutation header must set bulk=true");
        if (mutation.connection == null)
            throw new MutationValidationException("Mutation has no connection");
        provider = GrokConnect.providerManager.getByName(mutation.connection.dataSource);
        if (provider == null)
            throw new MutationValidationException("Unknown data source: " + mutation.connection.dataSource);
        if (!provider.descriptor.supportsBulkInsert)
            throw new UnsupportedOperationException("Bulk insert is not supported for provider " + provider.descriptor.type);
    }

    /** Opens the connection, starts the transaction and creates the loader. */
    public void start() throws SQLException, GrokConnectException {
        if (GrokConnectUtil.isNotEmpty(mutation.catalog))
            mutation.connection.parameters.put("db", mutation.catalog);
        try {
            connection = provider.getConnection(mutation.connection);
            provider.configureAutoCommit(connection);
            loader = provider.createBulkLoader(connection, mutation);
        } catch (SQLException | RuntimeException e) {
            abort();
            throw e;
        }
    }

    public void feed(byte[] csvChunk) throws Exception {
        if (finished || loader == null)
            throw new MutationValidationException("Mutation session is not active");
        loader.feed(csvChunk);
    }

    /** Flushes the loader, commits and closes the connection. */
    public MutationResult finish() throws Exception {
        if (finished || loader == null)
            throw new MutationValidationException("Mutation session is not active");
        MutationResult result;
        try {
            result = loader.finish();
            if (connection != null && !connection.getAutoCommit())
                connection.commit();
        } catch (Exception e) {
            abort();
            throw e;
        }
        finished = true;
        closeConnectionQuietly();
        return result;
    }

    /** Rolls back and releases everything; idempotent and never throws. */
    public void abort() {
        if (finished)
            return;
        finished = true;
        if (loader != null)
            try {
                loader.abort();
            } catch (Throwable t) {
                LOGGER.warn("Failed to abort bulk loader", t);
            }
        if (connection != null)
            provider.rollbackQuietly(connection);
        closeConnectionQuietly();
    }

    private void closeConnectionQuietly() {
        if (connection == null)
            return;
        try {
            if (!connection.isClosed())
                connection.close();
        } catch (SQLException e) {
            LOGGER.warn("Failed to close mutation connection", e);
        }
    }
}
