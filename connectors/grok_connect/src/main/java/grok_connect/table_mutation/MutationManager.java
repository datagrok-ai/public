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
 * lifecycle: {@link #start()} opens the connection and loader, {@link #feed} consumes decoded d42
 * chunks, {@link #finish()} flushes and commits, {@link #abort()} rolls back. The connection is always
 * closed before it returns to the pool (GROK-20323 contract). Unlike QueryManager this never
 * commits on close — commit happens only in {@link #finish()}.
 */
public class MutationManager {
    private static final Logger LOGGER = LoggerFactory.getLogger(MutationManager.class);
    private static final String CSV = "csv";
    private static final String D42 = "d42";

    public final JdbcDataProvider provider;
    private final InsertRows mutation;
    private final String mainCallId;
    private Connection connection;
    private BulkLoader loader;
    private boolean finished;
    private String createLeftoverNote; // §3.4.3 honesty: set once the create leg ran on a non-transactional-DDL dialect

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
        mainCallId = (String) call.aux.get("mainCallId");
        mutation = (InsertRows) call.func;
        if (!mutation.bulk)
            throw new MutationValidationException("Bulk mutation header must set bulk=true");
        if (mutation.connection == null)
            throw new MutationValidationException("Mutation has no connection");
        if (CSV.equals(mutation.payloadFormat))
            throw new MutationValidationException("CSV bulk payload is no longer supported; this GrokConnect requires "
                    + "payloadFormat='d42'. Datlas and GrokConnect must be upgraded together (java-d42-reader).");
        if (!D42.equals(mutation.payloadFormat))
            throw new MutationValidationException("Unknown payloadFormat '" + mutation.payloadFormat + "' (expected d42)");
        provider = GrokConnect.providerManager.getByName(mutation.connection.dataSource);
        if (provider == null)
            throw new MutationValidationException("Unknown data source: " + mutation.connection.dataSource);
        if (!provider.descriptor.supportsBulkInsert)
            throw new UnsupportedOperationException("Bulk insert is not supported for provider " + provider.descriptor.type);
        if ("create".equals(mutation.mode) && !provider.descriptor.supportsDdl)
            throw new UnsupportedOperationException("DDL is not supported for provider " + provider.descriptor.type);
    }

    /** Opens the connection, starts the transaction, runs the mode=create leg (§3.4.3) and creates the loader. */
    public void start() throws SQLException, GrokConnectException {
        if (GrokConnectUtil.isNotEmpty(mutation.catalog))
            mutation.connection.parameters.put("db", mutation.catalog);
        try {
            connection = provider.getConnection(mutation.connection);
            provider.configureAutoCommit(connection);
            if ("create".equals(mutation.mode)) {
                // mode=create (§3.4.3): create the derived all-nullable table BEFORE the loader prepares its
                // INSERT/COPY (and before MUTATION READY); post-creation the load is a plain insert, so the
                // provider's loader routing (the Postgres COPY guard) needs no create case.
                DdlRunner.execute(provider, connection, DdlRunner.deriveCreateTable(mutation), mainCallId);
                if (!provider.descriptor.supportsTransactionalDdl)
                    createLeftoverNote = DdlRunner.createLeftoverNote(provider, mutation);
                mutation.mode = "insert";
            }
            loader = provider.createBulkLoader(connection, mutation);
        } catch (SQLException | RuntimeException e) {
            abort();
            // loader creation may fail after the auto-committed create leg — the note must ride along
            Exception noted = withCreateLeftoverNote(e);
            if (noted instanceof SQLException)
                throw (SQLException) noted;
            throw (RuntimeException) noted;
        }
    }

    public void feed(byte[] bytes) throws Exception {
        if (finished || loader == null)
            throw new MutationValidationException("Mutation session is not active");
        try {
            DataFrame chunk = DataFrame.fromByteArray(bytes);
            validateChunkSchema(chunk);
            loader.feed(chunk);
        } catch (Exception e) {
            throw withCreateLeftoverNote(e);
        }
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
            if (rollback) {
                provider.rollbackQuietly(connection);
                if (createLeftoverNote != null) // the DML rolled back, but the auto-committed CREATE did not
                    result.errorMessage = (result.errorMessage == null
                            ? "Bulk load failed and was rolled back" : result.errorMessage) + createLeftoverNote;
            }
            else if (connection != null && !connection.getAutoCommit())
                connection.commit();
        } catch (Exception e) {
            abort();
            throw withCreateLeftoverNote(e);
        }
        finished = true;
        closeConnectionQuietly();
        return result;
    }

    /** §3.4.3 honesty: after the create leg on a non-transactional-DDL dialect, a load failure must say the table remains. */
    private Exception withCreateLeftoverNote(Exception e) {
        if (createLeftoverNote == null)
            return e;
        if (e instanceof MutationValidationException)
            return new MutationValidationException(e.getMessage() + createLeftoverNote,
                    ((MutationValidationException) e).refusals);
        if (e instanceof SQLException) {
            SQLException s = (SQLException) e;
            SQLException wrapped = new SQLException(s.getMessage() + createLeftoverNote, s.getSQLState(), s.getErrorCode(), s);
            wrapped.setNextException(s); // keep driver-specific chains (server-error walkers) intact
            return wrapped;
        }
        return new RuntimeException(e.getMessage() + createLeftoverNote, e);
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
