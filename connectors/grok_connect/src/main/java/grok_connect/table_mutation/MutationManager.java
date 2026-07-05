package grok_connect.table_mutation;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.List;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.GrokConnectUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.DataFrame;

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
    private static final String CSV = "csv";
    private static final String D42 = "d42";

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
        if (!CSV.equals(mutation.payloadFormat) && !D42.equals(mutation.payloadFormat))
            throw new MutationValidationException("Unknown payloadFormat '" + mutation.payloadFormat + "' (expected csv|d42)");
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

    public void feed(byte[] bytes) throws Exception {
        if (finished || loader == null)
            throw new MutationValidationException("Mutation session is not active");
        if (D42.equals(mutation.payloadFormat)) {
            DataFrame chunk = DataFrame.fromByteArray(bytes);
            validateChunkSchema(chunk);
            loader.feed(chunk);
        }
        else
            loader.feedCsv(bytes);
    }

    /**
     * Rejects a chunk whose columns do not match the header {@code columns}/{@code columnTypes} exactly
     * (names and dg types, order included) — an integrity check the CSV transport never had. The Datlas
     * sender controls the column order, so no order leniency. A mismatch throws
     * {@link MutationValidationException}, surfacing through the existing {@code ERROR:} + rollback path.
     */
    private void validateChunkSchema(DataFrame chunk) {
        List<String> columns = mutation.columns;
        List<String> types = mutation.columnTypes;
        if (columns == null || types == null || types.size() != columns.size())
            throw new MutationValidationException("Bulk mutation requires columnTypes parallel to columns");
        if (chunk.getColumnCount() != columns.size())
            throw new MutationValidationException("d42 chunk has " + chunk.getColumnCount()
                    + " column(s), expected " + columns.size());
        for (int i = 0; i < columns.size(); i++) {
            Column<?> col = chunk.getColumn(i);
            if (!columns.get(i).equals(col.getName()))
                throw new MutationValidationException("d42 chunk column " + i + " is '" + col.getName()
                        + "', expected '" + columns.get(i) + "'");
            if (!types.get(i).equals(col.getType()))
                throw new MutationValidationException("d42 chunk column '" + col.getName() + "' has type '"
                        + col.getType() + "', expected '" + types.get(i) + "'");
        }
    }

    /**
     * Flushes the loader, then commits — or, for an all-or-nothing load that hit a constraint error,
     * rolls the whole transaction back (connector-writes WO-6) — and closes the connection. Partial-mode
     * loads (allOrNothing=false) always commit the surviving rows; the per-row errors travel in the result.
     */
    public MutationResult finish() throws Exception {
        if (finished || loader == null)
            throw new MutationValidationException("Mutation session is not active");
        MutationResult result;
        try {
            result = loader.finish();
            boolean rollback = mutation.allOrNothing && result.errorCount != null && result.errorCount > 0;
            if (rollback)
                provider.rollbackQuietly(connection);
            else if (connection != null && !connection.getAutoCommit())
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
