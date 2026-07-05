package grok_connect.table_mutation;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import grok_connect.connectors_info.FuncCall;
import grok_connect.connectors_info.FuncParam;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.GrokConnectException;
import grok_connect.utils.GrokConnectUtil;
import grok_connect.utils.QueryCancelledByUser;
import grok_connect.utils.QueryMonitor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Types;

/**
 * Executes a TableMutation (or MutationBatch) on a JDBC provider inside one transaction:
 * commit on full success, rollback on any error, always before the connection returns to the
 * pool (GROK-20323 contract). SQL execution errors are reported inside the returned
 * {@link MutationResult} (errorMessage + errors), mirroring how /query reports them;
 * payload problems throw {@link MutationValidationException}.
 */
public class MutationRunner {
    private static final Logger LOGGER = LoggerFactory.getLogger(MutationRunner.class);
    private static final int BATCH_SIZE = 1000;
    /** Generic SQLState-independent error code; refined per-state by WO-6's SqlStateMapper. */
    private static final String DB_ERROR_CODE = "db-error";

    public static MutationResult execute(JdbcDataProvider provider, FuncCall call) throws GrokConnectException, QueryCancelledByUser {
        TableMutation mutation = (TableMutation) call.func;
        if (mutation.connection == null)
            throw new MutationValidationException("Mutation has no connection");
        if (GrokConnectUtil.isNotEmpty(mutation.catalog))
            mutation.connection.parameters.put("db", mutation.catalog);
        String mainCallId = (String) call.aux.get("mainCallId");
        MutationResult result = new MutationResult();
        List<TableMutation> operations = mutation instanceof MutationBatch
                ? ((MutationBatch) mutation).operations : Collections.singletonList(mutation);
        if (operations == null || operations.isEmpty())
            throw new MutationValidationException("MutationBatch requires a non-empty operations list");
        if (mutation instanceof MutationBatch)
            result.perStatement = new ArrayList<>();
        Connection connection = null;
        int currentStatement = 0;
        try {
            connection = provider.getConnection(mutation.connection);
            provider.configureAutoCommit(connection);
            for (int i = 0; i < operations.size(); i++) {
                currentStatement = i;
                int affected = executeOperation(provider, connection, operations.get(i), mainCallId);
                result.affectedRows += affected;
                if (result.perStatement != null) {
                    PerStatementResult statementResult = new PerStatementResult();
                    statementResult.statementIndex = i;
                    statementResult.affectedRows = affected;
                    result.perStatement.add(statementResult);
                }
            }
            if (!connection.getAutoCommit())
                connection.commit();
            return result;
        } catch (SQLException e) {
            provider.rollbackQuietly(connection);
            if (QueryMonitor.getInstance().checkCancelledId(mainCallId))
                throw new QueryCancelledByUser();
            LOGGER.info("Mutation failed and was rolled back", e);
            result.errorMessage = e.getMessage();
            RowError error = new RowError();
            error.index = currentStatement;
            error.code = DB_ERROR_CODE;
            error.message = e.getMessage();
            result.errors = Collections.singletonList(error);
            result.errorCount = 1;
            return result;
        } catch (RuntimeException e) {
            provider.rollbackQuietly(connection);
            throw e;
        } finally {
            if (connection != null)
                try {
                    connection.close();
                } catch (SQLException e) {
                    LOGGER.warn("Failed to close connection", e);
                }
        }
    }

    private static int executeOperation(JdbcDataProvider provider, Connection connection, TableMutation m, String mainCallId) throws SQLException {
        if (m instanceof MutationBatch) {
            List<TableMutation> operations = ((MutationBatch) m).operations;
            if (operations == null || operations.isEmpty())
                throw new MutationValidationException("MutationBatch requires a non-empty operations list");
            int affected = 0;
            for (TableMutation operation : operations)
                affected += executeOperation(provider, connection, operation, mainCallId);
            return affected;
        }
        if (m instanceof UpsertRows) // WO-4 adds upsertSql; until then upsert is a capability error
            throw new UnsupportedOperationException("Upsert is not supported yet for provider " + provider.descriptor.type);
        if (m instanceof InsertRows)
            return executeInsert(provider, connection, (InsertRows) m, mainCallId);
        if (m instanceof UpdateRows)
            return executeUpdate(provider, connection, (UpdateRows) m, mainCallId);
        if (m instanceof DeleteRows)
            return executeDelete(provider, connection, (DeleteRows) m, mainCallId);
        throw new MutationValidationException("Unsupported mutation type: " + m.type);
    }

    private static int executeInsert(JdbcDataProvider provider, Connection connection, InsertRows m, String mainCallId) throws SQLException {
        if (GrokConnectUtil.isEmpty(m.mode)) // hand-built payloads may omit it; the contract default is insert
            m.mode = "insert";
        if (m.bulk && m.rows == null) // WO-5 adds the streamed payload
            throw new UnsupportedOperationException("Bulk insert streaming is not supported yet");
        String sql = provider.insertSql(m);
        if (m.columnTypes == null || m.columnTypes.size() != m.columns.size())
            throw new MutationValidationException("InsertRows requires columnTypes parallel to columns");
        if (m.rows == null || m.rows.isEmpty())
            throw new MutationValidationException("InsertRows requires a non-empty rows list");
        LOGGER.info("Mutation before execution: {}", sql);
        QueryMonitor queryMonitor = QueryMonitor.getInstance();
        try (PreparedStatement statement = connection.prepareStatement(sql)) {
            queryMonitor.addNewStatement(mainCallId, statement);
            try {
                int affected = 0;
                int pending = 0;
                for (int r = 0; r < m.rows.size(); r++) {
                    List<Object> row = m.rows.get(r);
                    if (row == null || row.size() != m.columns.size())
                        throw new MutationValidationException("Row " + r + " size does not match the columns list");
                    for (int c = 0; c < row.size(); c++)
                        bindValue(provider, statement, c + 1, row.get(c), m.columnTypes.get(c));
                    statement.addBatch();
                    if (++pending == BATCH_SIZE) {
                        affected += sumBatchCounts(statement.executeBatch());
                        pending = 0;
                    }
                }
                if (pending > 0)
                    affected += sumBatchCounts(statement.executeBatch());
                return affected;
            } finally {
                queryMonitor.removeStatement(mainCallId);
            }
        }
    }

    private static int executeUpdate(JdbcDataProvider provider, Connection connection, UpdateRows m, String mainCallId) throws SQLException {
        List<FuncParam> whereParams = new ArrayList<>();
        String sql = provider.updateSql(m, whereParams);
        return executeUpdateOrDelete(provider, connection, m, sql, whereParams, mainCallId);
    }

    private static int executeDelete(JdbcDataProvider provider, Connection connection, DeleteRows m, String mainCallId) throws SQLException {
        List<FuncParam> whereParams = new ArrayList<>();
        String sql = provider.deleteSql(m, whereParams);
        return executeUpdateOrDelete(provider, connection, m, sql, whereParams, mainCallId);
    }

    private static int executeUpdateOrDelete(JdbcDataProvider provider, Connection connection, TableMutation m,
                                             String sql, List<FuncParam> whereParams, String mainCallId) throws SQLException {
        // getParameterNames resolves the @-named predicate params against m.params
        m.params = whereParams;
        StringBuilder queryBuffer = new StringBuilder();
        List<String> names = provider.getParameterNames(sql, m, queryBuffer);
        sql = queryBuffer.toString();
        if (names.size() != whereParams.size())
            throw new MutationValidationException("Predicate parameter mismatch: " + whereParams.size()
                    + " collected, " + names.size() + " referenced in SQL");
        LOGGER.info("Mutation before execution: {}", sql);
        QueryMonitor queryMonitor = QueryMonitor.getInstance();
        try (PreparedStatement statement = connection.prepareStatement(sql)) {
            queryMonitor.addNewStatement(mainCallId, statement);
            try {
                int idx = 1;
                if (m instanceof UpdateRows) {
                    UpdateRows update = (UpdateRows) m;
                    for (int i = 0; i < update.setColumns.size(); i++)
                        bindValue(provider, statement, idx++, update.setValues.get(i), update.setTypes.get(i));
                }
                // pattern converters emit @-placeholders in the same order they collect params,
                // so binding is positional — immune to duplicate param names across predicates
                for (FuncParam param : whereParams)
                    bindValue(provider, statement, idx++, param.value, param.propertyType);
                return statement.executeUpdate();
            } finally {
                queryMonitor.removeStatement(mainCallId);
            }
        }
    }

    /** Binds a value coerced by its declared dg type (JSON numbers arrive as Double, bigints as strings). */
    public static void bindValue(JdbcDataProvider provider, PreparedStatement statement, int idx, Object value, String dgType) throws SQLException {
        if (dgType == null)
            throw new MutationValidationException("Missing dg type for parameter " + idx);
        switch (dgType) {
            case Types.INT:
                if (value == null)
                    statement.setNull(idx, java.sql.Types.INTEGER);
                else
                    statement.setInt(idx, value instanceof Number ? ((Number) value).intValue() : Integer.parseInt(value.toString()));
                break;
            case Types.BIG_INT:
                if (value == null)
                    statement.setNull(idx, java.sql.Types.BIGINT);
                else
                    statement.setLong(idx, value instanceof Number ? ((Number) value).longValue() : Long.parseLong(value.toString()));
                break;
            case Types.NUM:
            case Types.FLOAT:
                if (value == null)
                    statement.setNull(idx, java.sql.Types.DOUBLE);
                else
                    statement.setDouble(idx, value instanceof Number ? ((Number) value).doubleValue() : Double.parseDouble(value.toString()));
                break;
            case Types.BOOL:
                if (value == null)
                    statement.setNull(idx, java.sql.Types.BOOLEAN);
                else
                    statement.setBoolean(idx, value instanceof Boolean ? (Boolean) value : Boolean.parseBoolean(value.toString()));
                break;
            case Types.DATE_TIME:
                provider.setDateTimeValue(new FuncParam(Types.DATE_TIME, "", value == null ? null : value.toString()), statement, idx);
                break;
            case Types.STRING:
                if (value == null)
                    statement.setNull(idx, java.sql.Types.VARCHAR);
                else
                    statement.setString(idx, value.toString());
                break;
            default:
                throw new MutationValidationException("Unsupported dg type '" + dgType + "' for parameter " + idx);
        }
    }

    private static int sumBatchCounts(int[] counts) {
        int sum = 0;
        for (int count : counts)
            sum += count == Statement.SUCCESS_NO_INFO ? 1 : Math.max(count, 0);
        return sum;
    }
}
