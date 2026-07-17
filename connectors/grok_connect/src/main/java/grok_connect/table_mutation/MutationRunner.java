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
        if (mutation instanceof MutationBatch) {
            refuseSingleOpMembers(operations); // DDL and mode=create execute one operation per request (§3.4.1)
            result.perStatement = new ArrayList<>();
        }
        boolean createMode = mutation instanceof InsertRows && "create".equals(((InsertRows) mutation).mode);
        // dryRun passes the SAME capability gates as a real run (§3.4.2), so the gate precedes the dryRun branch
        if ((mutation instanceof DdlMutation || createMode) && !provider.descriptor.supportsDdl)
            throw new UnsupportedOperationException("DDL is not supported for provider " + provider.descriptor.type);
        if (mutation.dryRun && !(mutation instanceof DdlMutation)) {
            // DML dryRun: statement emission only — no binding, no pre-checks, no connection (§3.4.2);
            // create mode travels inline here too (bulk=false, rows=null — row validation is skipped, §3.4.3)
            result.plan = new MutationPlan();
            result.plan.statements = new ArrayList<>();
            for (TableMutation operation : operations)
                collectDmlStatements(provider, operation, result.plan.statements);
            if (createMode)
                result.plan.transactionalDdl = provider.descriptor.supportsTransactionalDdl;
            return result;
        }
        Connection connection = null;
        int currentStatement = 0;
        try {
            connection = provider.getConnection(mutation.connection);
            provider.configureAutoCommit(connection);
            if (mutation instanceof DdlMutation) {
                DdlMutation ddl = (DdlMutation) mutation;
                if (mutation.dryRun) {
                    result.plan = DdlRunner.plan(provider, connection, ddl);
                    provider.rollbackQuietly(connection); // pre-check SELECTs joined the transaction; nothing ran
                    return result;
                }
                result.plan = DdlRunner.execute(provider, connection, ddl, mainCallId); // affectedRows stays 0
            }
            else
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
            result.errors = Collections.singletonList(SqlStateMapper.toRowError(provider, currentStatement, e));
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
        if (m instanceof UpsertRows) // must precede InsertRows: UpsertRows extends InsertRows
            return executeUpsert(provider, connection, (UpsertRows) m, mainCallId);
        if (m instanceof InsertRows)
            return executeInsert(provider, connection, (InsertRows) m, mainCallId);
        if (m instanceof UpdateRows)
            return executeUpdate(provider, connection, (UpdateRows) m, mainCallId);
        if (m instanceof DeleteRows)
            return executeDelete(provider, connection, (DeleteRows) m, mainCallId);
        throw new MutationValidationException("Unsupported mutation type: " + m.type);
    }

    /**
     * DDL is strictly single-op-per-request (ARCHITECTURE §3.4.1, resolved open question §12.1):
     * many DBs auto-commit DDL, so batch atomicity would be a lie. Datlas enforces the same rule
     * fail-fast; this is the defense in depth.
     */
    private static void refuseSingleOpMembers(List<TableMutation> operations) {
        for (TableMutation operation : operations) {
            if (operation instanceof DdlMutation
                    || (operation instanceof InsertRows && "create".equals(((InsertRows) operation).mode)))
                throw new MutationValidationException("DDL executes one operation per request");
            if (operation instanceof MutationBatch && ((MutationBatch) operation).operations != null)
                refuseSingleOpMembers(((MutationBatch) operation).operations);
        }
    }

    /**
     * DML dryRun statements (§3.4.2): the same emitters the execute path uses, with predicate
     * params compiled to placeholders — no binding, no pre-checks. Upserts are emitted in the
     * one-row shape (chunked dialects repeat the value tuple at execution time).
     */
    private static void collectDmlStatements(JdbcDataProvider provider, TableMutation m, List<String> statements) {
        if (m instanceof MutationBatch) {
            for (TableMutation operation : ((MutationBatch) m).operations)
                collectDmlStatements(provider, operation, statements);
            return;
        }
        if (m instanceof UpdateRows || m instanceof DeleteRows) {
            List<FuncParam> whereParams = new ArrayList<>();
            String sql = m instanceof UpdateRows
                    ? provider.updateSql((UpdateRows) m, whereParams)
                    : provider.deleteSql((DeleteRows) m, whereParams);
            m.params = whereParams;
            StringBuilder queryBuffer = new StringBuilder();
            provider.getParameterNames(sql, m, queryBuffer);
            statements.add(queryBuffer.toString());
            return;
        }
        if (m instanceof UpsertRows) {
            requireUpsertSupport(provider);
            statements.add(provider.upsertSql((UpsertRows) m, 1));
            return;
        }
        if (m instanceof InsertRows) {
            InsertRows insert = (InsertRows) m;
            if (GrokConnectUtil.isEmpty(insert.mode))
                insert.mode = "insert";
            switch (insert.mode) {
                case "insert":
                    statements.add(provider.insertSql(insert));
                    return;
                case "upsert":
                    requireUpsertSupport(provider);
                    statements.add(provider.upsertSql(new UpsertRows(insert), 1));
                    return;
                case "update":
                    statements.add(provider.updateByKeySql(insert));
                    return;
                case "create":
                    // the CREATE statements + the insert shape (§3.4.3) — same emitters as an explicit CreateTable
                    statements.addAll(provider.createTableSql(DdlRunner.deriveCreateTable(insert)));
                    statements.add(provider.insertSql(insert));
                    return;
                default:
                    throw new MutationValidationException("Unsupported insert mode: " + insert.mode);
            }
        }
        throw new MutationValidationException("Unsupported mutation type: " + m.type);
    }

    /** dryRun passes the same capability gates as a real run (§3.4.2) — the descriptor flag, not the emitter, is the authority. */
    private static void requireUpsertSupport(JdbcDataProvider provider) {
        if (!provider.descriptor.supportsUpsert)
            throw new UnsupportedOperationException("Upsert is not supported for provider " + provider.descriptor.type);
    }

    private static int executeInsert(JdbcDataProvider provider, Connection connection, InsertRows m, String mainCallId) throws SQLException {
        if (GrokConnectUtil.isEmpty(m.mode)) // hand-built payloads may omit it; the contract default is insert
            m.mode = "insert";
        // an InsertRows carrying mode=upsert/update/create routes to the matching engine (bulk parity, WO-6);
        // any other mode fails loud rather than silently doing a plain insert
        switch (m.mode) {
            case "insert":
                return executePlainInsert(provider, connection, m, mainCallId);
            case "upsert":
                return executeUpsert(provider, connection, m instanceof UpsertRows ? (UpsertRows) m : new UpsertRows(m), mainCallId);
            case "update":
                return executeUpdateByKey(provider, connection, m, mainCallId);
            case "create":
                return executeCreate(provider, connection, m, mainCallId);
            default:
                throw new MutationValidationException("Unsupported insert mode: " + m.mode);
        }
    }

    /**
     * mode=create (ARCHITECTURE §3.4.3): derives an all-nullable key-less CreateTable from
     * columns/columnTypes, creates it on the same connection through {@link DdlRunner} (the single
     * emission path — the same emitters an explicit CreateTable uses), then runs the plain insert
     * engine. Payload validation runs BEFORE the CREATE so a bad payload never leaves a table behind.
     * On transactional-DDL dialects a failed load rolls the CREATE back with the transaction; elsewhere
     * the created table remains, empty — reported honestly in the error message.
     */
    private static int executeCreate(JdbcDataProvider provider, Connection connection, InsertRows m, String mainCallId) throws SQLException {
        CreateTable create = DdlRunner.deriveCreateTable(m); // validates columns/columnTypes before any DDL
        if (m.bulk && m.rows == null) // the streamed create path runs through MutationManager, not /mutate
            throw new UnsupportedOperationException("Bulk insert streaming is not supported yet");
        if (m.rows == null || m.rows.isEmpty())
            throw new MutationValidationException("InsertRows requires a non-empty rows list");
        DdlRunner.execute(provider, connection, create, mainCallId); // additive — never destructive, no confirm
        try {
            return executePlainInsert(provider, connection, m, mainCallId);
        } catch (SQLException e) {
            throw withCreateLeftoverNote(provider, m, e);
        } catch (MutationValidationException e) {
            if (provider.descriptor.supportsTransactionalDdl)
                throw e;
            throw new MutationValidationException(e.getMessage() + DdlRunner.createLeftoverNote(provider, m), e.refusals);
        }
    }

    /** §3.4.3 honesty: on a non-transactional-DDL dialect the CREATE has already committed — a failed load leaves the empty table behind, and the error must say so. */
    private static SQLException withCreateLeftoverNote(JdbcDataProvider provider, InsertRows m, SQLException e) {
        if (provider.descriptor.supportsTransactionalDdl)
            return e;
        SQLException wrapped = new SQLException(e.getMessage() + DdlRunner.createLeftoverNote(provider, m),
                e.getSQLState(), e.getErrorCode(), e);
        wrapped.setNextException(e); // keep driver-specific chains (server-error walkers) intact
        return wrapped;
    }

    private static int executePlainInsert(JdbcDataProvider provider, Connection connection, InsertRows m, String mainCallId) throws SQLException {
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

    private static int executeUpsert(JdbcDataProvider provider, Connection connection, UpsertRows m, String mainCallId) throws SQLException {
        requireUpsertSupport(provider);
        if (m.bulk && m.rows == null) // WO-5 adds the streamed payload
            throw new UnsupportedOperationException("Bulk upsert streaming is not supported yet");
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("UpsertRows requires a non-empty columns list");
        if (m.columnTypes == null || m.columnTypes.size() != m.columns.size())
            throw new MutationValidationException("UpsertRows requires columnTypes parallel to columns");
        if (m.matchKeys == null || m.matchKeys.isEmpty())
            throw new MutationValidationException("UpsertRows requires a non-empty matchKeys list");
        if (m.rows == null || m.rows.isEmpty())
            throw new MutationValidationException("UpsertRows requires a non-empty rows list");
        int chunkRows = provider.upsertBatchRows(m.columns.size());
        return chunkRows <= 1
                ? executeUpsertBatched(provider, connection, m, mainCallId)
                : executeUpsertChunked(provider, connection, m, chunkRows, mainCallId);
    }

    /**
     * Key-based batched UPDATE for inline {@code mode == "update"} rows: one row of placeholders per row.
     * Unlike the bulk loader, the inline path does not report missing keys in phase A — a row whose key
     * matches nothing (update count 0) is silently a no-op, not a {@code missing} RowError.
     */
    private static int executeUpdateByKey(JdbcDataProvider provider, Connection connection, InsertRows m, String mainCallId) throws SQLException {
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("Update mode requires a non-empty columns list");
        if (m.columnTypes == null || m.columnTypes.size() != m.columns.size())
            throw new MutationValidationException("Update mode requires columnTypes parallel to columns");
        if (m.rows == null || m.rows.isEmpty())
            throw new MutationValidationException("Update mode requires a non-empty rows list");
        String sql = provider.updateByKeySql(m);
        int[] bindOrder = provider.updateByKeyBindOrder(m);
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
                    for (int i = 0; i < bindOrder.length; i++)
                        bindValue(provider, statement, i + 1, row.get(bindOrder[i]), m.columnTypes.get(bindOrder[i]));
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

    /** addBatch single-row upsert (Postgres/MySQL: ON CONFLICT / ON DUPLICATE KEY; Oracle: MERGE FROM dual). */
    private static int executeUpsertBatched(JdbcDataProvider provider, Connection connection, UpsertRows m, String mainCallId) throws SQLException {
        String sql = provider.upsertSql(m, 1);
        LOGGER.info("Mutation before execution: {}", sql);
        QueryMonitor queryMonitor = QueryMonitor.getInstance();
        try (PreparedStatement statement = connection.prepareStatement(sql)) {
            queryMonitor.addNewStatement(mainCallId, statement);
            try {
                int affected = 0;
                int pending = 0;
                for (int r = 0; r < m.rows.size(); r++) {
                    bindRow(provider, statement, m, r, 1);
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

    /** Multi-row VALUES chunks, one executeUpdate per chunk (MS SQL / Snowflake MERGE-over-VALUES). */
    private static int executeUpsertChunked(JdbcDataProvider provider, Connection connection, UpsertRows m, int chunkRows, String mainCallId) throws SQLException {
        QueryMonitor queryMonitor = QueryMonitor.getInstance();
        int total = m.rows.size();
        int affected = 0;
        for (int start = 0; start < total; start += chunkRows) {
            int end = Math.min(start + chunkRows, total);
            String sql = provider.upsertSql(m, end - start);
            LOGGER.info("Mutation before execution: {}", sql);
            try (PreparedStatement statement = connection.prepareStatement(sql)) {
                queryMonitor.addNewStatement(mainCallId, statement);
                try {
                    int idx = 1;
                    for (int r = start; r < end; r++)
                        idx = bindRow(provider, statement, m, r, idx);
                    affected += Math.max(statement.executeUpdate(), 0);
                } finally {
                    queryMonitor.removeStatement(mainCallId);
                }
            }
        }
        return affected;
    }

    /** Binds row {@code r} starting at parameter {@code idx}; returns the next free parameter index. */
    private static int bindRow(JdbcDataProvider provider, PreparedStatement statement, InsertRows m, int r, int idx) throws SQLException {
        List<Object> row = m.rows.get(r);
        if (row == null || row.size() != m.columns.size())
            throw new MutationValidationException("Row " + r + " size does not match the columns list");
        for (int c = 0; c < row.size(); c++)
            bindValue(provider, statement, idx++, row.get(c), m.columnTypes.get(c));
        return idx;
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

    /**
     * Binds cell {@code row} of a decoded d42 {@code col} by its declared dg type (java-d42-reader WO-4):
     * {@code isNone} → {@code setNull} with the matching {@link java.sql.Types}; datetime epoch-µs →
     * {@link java.sql.Timestamp} with a UTC calendar (no string round-trip); ints/bigints/doubles/bools/
     * strings via the typed getters. The delegating {@link #bindValue} stays for inline {@code /mutate} rows.
     */
    public static void bindColumnValue(JdbcDataProvider provider, PreparedStatement statement, int idx,
                                       serialization.Column<?> col, int row, String dgType) throws SQLException {
        if (dgType == null)
            throw new MutationValidationException("Missing dg type for parameter " + idx);
        if (col.isNone(row)) {
            statement.setNull(idx, sqlTypeFor(dgType, idx));
            return;
        }
        switch (dgType) {
            case Types.INT:
                statement.setInt(idx, ((serialization.IntColumn) col).get(row));
                break;
            case Types.BIG_INT:
                statement.setLong(idx, Long.parseLong(col.get(row).toString()));
                break;
            case Types.NUM:
            case Types.FLOAT:
                statement.setDouble(idx, ((serialization.FloatColumn) col).getDouble(row));
                break;
            case Types.BOOL:
                statement.setBoolean(idx, ((serialization.BoolColumn) col).get(row));
                break;
            case Types.DATE_TIME:
                statement.setTimestamp(idx, microsToTimestamp(((serialization.DateTimeColumn) col).get(row)),
                        java.util.Calendar.getInstance(java.util.TimeZone.getTimeZone("UTC")));
                break;
            case Types.STRING:
                statement.setString(idx, col.get(row).toString());
                break;
            default:
                throw new MutationValidationException("Unsupported dg type '" + dgType + "' for parameter " + idx);
        }
    }

    private static int sqlTypeFor(String dgType, int idx) {
        switch (dgType) {
            case Types.INT: return java.sql.Types.INTEGER;
            case Types.BIG_INT: return java.sql.Types.BIGINT;
            case Types.NUM:
            case Types.FLOAT: return java.sql.Types.DOUBLE;
            case Types.BOOL: return java.sql.Types.BOOLEAN;
            case Types.DATE_TIME: return java.sql.Types.TIMESTAMP;
            case Types.STRING: return java.sql.Types.VARCHAR;
            default:
                throw new MutationValidationException("Unsupported dg type '" + dgType + "' for parameter " + idx);
        }
    }

    /** Epoch-µs (UTC) → Timestamp: whole-second base + µs-of-second in nanos (floor math, negative-safe). */
    static java.sql.Timestamp microsToTimestamp(double micros) {
        long us = (long) micros;
        long seconds = Math.floorDiv(us, 1_000_000L);
        long microOfSecond = Math.floorMod(us, 1_000_000L);
        java.sql.Timestamp ts = new java.sql.Timestamp(seconds * 1000L);
        ts.setNanos((int) (microOfSecond * 1000L));
        return ts;
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

    static int sumBatchCounts(int[] counts) {
        int sum = 0;
        for (int count : counts)
            sum += count == Statement.SUCCESS_NO_INFO ? 1 : Math.max(count, 0);
        return sum;
    }
}
