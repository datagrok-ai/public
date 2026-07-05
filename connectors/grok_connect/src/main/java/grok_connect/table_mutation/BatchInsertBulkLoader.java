package grok_connect.table_mutation;

import java.sql.BatchUpdateException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Savepoint;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.GrokConnectUtil;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.DataFrame;

/**
 * Default bulk loader for any prepared-statement JDBC provider: binds each decoded d42 row via
 * {@link MutationRunner#bindColumnValue} using the declared dg types (the legacy CSV path,
 * {@link #feedCsv}, still parses with {@link CsvChunkParser} and binds via
 * {@link MutationRunner#bindValue}), and {@code executeBatch}es every {@value #BATCH_SIZE} rows and at
 * {@link #finish()}.
 *
 * <p>Batch semantics (connector-writes WO-6): the {@code mode} selects the emitted SQL —
 * {@code insert} (optionally {@code ON CONFLICT DO NOTHING} when {@code errorOnDuplicate == false}),
 * {@code upsert} (dialect upsert SQL), or {@code update} (key-based UPDATE); any other mode fails loud.
 * {@code allOrNothing == true} (default) runs each chunk atomically — the first constraint error stops
 * the load and the manager rolls the whole transaction back, reporting the exact failing row. With
 * {@code allOrNothing == false} each chunk is wrapped in a savepoint; on a batch error the chunk rolls
 * back to the savepoint and replays row-by-row (each under its own savepoint), keeping the survivors and
 * collecting a per-row {@link RowError} (capped at {@value #MAX_ERRORS}; {@code errorCount} stays exact).
 * Holds at most one pending chunk plus the parser's partial record, so heap stays flat.
 */
public class BatchInsertBulkLoader implements BulkLoader {
    private static final Logger LOGGER = LoggerFactory.getLogger(BatchInsertBulkLoader.class);
    private static final int BATCH_SIZE = 10_000;
    private static final int MAX_ERRORS = 1000;

    // interpretation of a per-row update count of 0, by mode
    private static final int ZERO_IGNORE = 0;   // not expected (plain insert / upsert)
    private static final int ZERO_SKIP = 1;     // insert ON CONFLICT DO NOTHING: a duplicate was skipped
    private static final int ZERO_MISSING = 2;  // update by key: no row matched the key columns

    private final JdbcDataProvider provider;
    private final Connection connection;
    private final InsertRows mutation;
    private final String mode;
    private final boolean allOrNothing;
    private final int columnCount;
    private final int[] bindOrder;
    private final int zeroPolicy;
    private final boolean savepoints;
    private final PreparedStatement statement;
    private final CsvChunkParser parser = new CsvChunkParser();
    // Rows added since the last flush, retained so a failed chunk can re-bind row-by-row on replay.
    // A d42 row captures its (chunk, row); a CSV row captures the parsed String list — the flush /
    // savepoint / replay engine below is identical for both.
    private final List<PendingRow> batch = new ArrayList<>();
    private final List<RowError> errors = new ArrayList<>();

    private int nextIndex;   // global index of the next row to be added
    private int affected;    // rows the DB actually changed
    private int inserted;
    private int updated;
    private int skipped;
    private int errorCount;
    private boolean stopped;  // an all-or-nothing failure occurred; ignore the remaining stream

    public BatchInsertBulkLoader(JdbcDataProvider provider, Connection connection, InsertRows m) throws SQLException {
        if (m.columns == null || m.columns.isEmpty())
            throw new MutationValidationException("Bulk insert requires a non-empty columns list");
        if (m.columnTypes == null || m.columnTypes.size() != m.columns.size())
            throw new MutationValidationException("Bulk insert requires columnTypes parallel to columns");
        this.provider = provider;
        this.connection = connection;
        this.mutation = m;
        this.columnCount = m.columns.size();
        this.mode = GrokConnectUtil.isEmpty(m.mode) ? "insert" : m.mode;
        this.allOrNothing = m.allOrNothing;
        String sql;
        switch (this.mode) {
            case "insert":
                sql = m.errorOnDuplicate ? provider.insertSql(m) : provider.insertIgnoreDuplicatesSql(m);
                bindOrder = identity(columnCount);
                zeroPolicy = m.errorOnDuplicate ? ZERO_IGNORE : ZERO_SKIP;
                break;
            case "upsert":
                UpsertRows upsert = m instanceof UpsertRows ? (UpsertRows) m : new UpsertRows(m);
                sql = provider.upsertSql(upsert, 1);
                bindOrder = identity(columnCount);
                zeroPolicy = ZERO_IGNORE;
                break;
            case "update":
                sql = provider.updateByKeySql(m);
                bindOrder = provider.updateByKeyBindOrder(m);
                zeroPolicy = ZERO_MISSING;
                break;
            default:
                throw new MutationValidationException("Unsupported bulk mutation mode: " + this.mode);
        }
        this.savepoints = connection.getMetaData().supportsSavepoints();
        if (!allOrNothing && !savepoints)
            throw new UnsupportedOperationException("Partial mode (allOrNothing=false) is not supported by provider "
                    + provider.descriptor.type);
        this.statement = connection.prepareStatement(sql);
    }

    private static int[] identity(int n) {
        int[] order = new int[n];
        for (int i = 0; i < n; i++)
            order[i] = i;
        return order;
    }

    @Override
    public void feed(DataFrame chunk) throws SQLException {
        if (chunk.getColumnCount() != columnCount)
            throw new MutationValidationException("d42 chunk has " + chunk.getColumnCount() + " column(s), expected " + columnCount);
        for (int r = 0; r < chunk.rowCount; r++) {
            final int row = r;
            addRow(() -> bindColumns(chunk, row));
        }
    }

    @Override
    public void feedCsv(byte[] csvChunk) throws SQLException {
        for (List<String> row : parser.feed(csvChunk))
            addCsvRow(row);
    }

    @Override
    public MutationResult finish() throws SQLException {
        for (List<String> row : parser.finish())  // no-op tail for the d42 path (the parser was never fed)
            addCsvRow(row);
        flush();
        closeQuietly();
        return buildResult();
    }

    @Override
    public void abort() {
        closeQuietly();
    }

    /** A row whose columns can be (re-)bound to {@link #statement} — from a d42 chunk or a parsed CSV record. */
    @FunctionalInterface
    private interface PendingRow {
        void bind() throws SQLException;
    }

    private void addCsvRow(List<String> row) throws SQLException {
        if (row.size() != columnCount)
            throw new MutationValidationException("CSV row has " + row.size() + " field(s), expected " + columnCount);
        addRow(() -> bindStrings(row));
    }

    private void addRow(PendingRow row) throws SQLException {
        nextIndex++;
        if (stopped)  // an all-or-nothing failure already poisoned the transaction; the manager rolls it back
            return;
        row.bind();
        statement.addBatch();
        batch.add(row);
        if (batch.size() == BATCH_SIZE)
            flush();
    }

    private void bindColumns(DataFrame chunk, int row) throws SQLException {
        for (int i = 0; i < bindOrder.length; i++) {
            int col = bindOrder[i];
            MutationRunner.bindColumnValue(provider, statement, i + 1, chunk.getColumn(col), row, mutation.columnTypes.get(col));
        }
    }

    private void bindStrings(List<String> row) throws SQLException {
        for (int i = 0; i < bindOrder.length; i++) {
            int col = bindOrder[i];
            MutationRunner.bindValue(provider, statement, i + 1, row.get(col), mutation.columnTypes.get(col));
        }
    }

    private void flush() throws SQLException {
        if (batch.isEmpty())
            return;
        int chunkStart = nextIndex - batch.size();
        if (allOrNothing)
            flushAtomic(chunkStart);
        else
            flushPartial(chunkStart);
        batch.clear();
    }

    private void flushAtomic(int chunkStart) throws SQLException {
        // Wrap the chunk in a savepoint (when supported) purely so a failure can be diagnosed: pgjdbc
        // reports the whole failed batch as EXECUTE_FAILED, so getUpdateCounts() cannot pinpoint the row.
        // On failure we roll back to the savepoint and replay to the first failing row (exact index), then
        // stop — the manager rolls the whole transaction back anyway.
        Savepoint savepoint = savepoints ? connection.setSavepoint() : null;
        try {
            applyCounts(statement.executeBatch(), chunkStart);
            if (savepoint != null)
                releaseQuietly(savepoint);
        } catch (BatchUpdateException e) {
            if (savepoint != null) {
                connection.rollback(savepoint);
                releaseQuietly(savepoint);
                statement.clearBatch();
                recordFirstFailure(chunkStart, e);
            }
            else
                recordFailingRowFromCounts(chunkStart, e);
            stopped = true; // the whole transaction is poisoned; the manager rolls it back
        }
    }

    private void flushPartial(int chunkStart) throws SQLException {
        Savepoint savepoint = connection.setSavepoint();
        try {
            applyCounts(statement.executeBatch(), chunkStart);
            releaseQuietly(savepoint);
        } catch (BatchUpdateException e) {
            connection.rollback(savepoint); // undo the whole failed chunk, then replay it row-by-row
            releaseQuietly(savepoint);
            statement.clearBatch();
            replayRowByRow(chunkStart);
        }
    }

    /** Re-runs each buffered row of a failed chunk under its own savepoint, keeping the survivors. */
    private void replayRowByRow(int chunkStart) throws SQLException {
        for (int i = 0; i < batch.size(); i++) {
            int rowIndex = chunkStart + i;
            Savepoint savepoint = connection.setSavepoint();
            try {
                batch.get(i).bind();
                classifyOne(statement.executeUpdate(), rowIndex);
                releaseQuietly(savepoint);
            } catch (SQLException e) {
                connection.rollback(savepoint);
                releaseQuietly(savepoint);
                addError(rowIndex, e);
            }
        }
    }

    /** Classifies a successful chunk's per-row update counts (array length varies by driver — read defensively). */
    private void applyCounts(int[] counts, int chunkStart) {
        for (int i = 0; i < batch.size(); i++) {
            int count = i < counts.length ? counts[i] : Statement.SUCCESS_NO_INFO;
            classifyOne(count, chunkStart + i);
        }
    }

    private void classifyOne(int count, int rowIndex) {
        if (count == 0) {
            if (zeroPolicy == ZERO_SKIP)
                skipped++;
            else if (zeroPolicy == ZERO_MISSING)
                addMissing(rowIndex);
            return;
        }
        int rows = count == Statement.SUCCESS_NO_INFO ? 1 : Math.max(count, 0);
        affected += rows;
        if (mode.equals("update"))
            updated += rows;
        else if (mode.equals("insert"))
            inserted += rows;
    }

    /**
     * Replays a rolled-back atomic chunk row-by-row to pinpoint the first failing row (exact index). If the
     * replay reproduces no error — the original failure was nondeterministic (deadlock 40P01 /
     * serialization 40001 that does not recur) — a statement-level error is still recorded from the
     * original {@code batchError} so {@code errorCount > 0} forces the manager's full rollback: never a
     * silent partial commit + false success.
     */
    private void recordFirstFailure(int chunkStart, BatchUpdateException batchError) throws SQLException {
        for (int i = 0; i < batch.size(); i++) {
            Savepoint savepoint = connection.setSavepoint();
            try {
                batch.get(i).bind();
                statement.executeUpdate();
                releaseQuietly(savepoint);
            } catch (SQLException e) {
                connection.rollback(savepoint);
                releaseQuietly(savepoint);
                addError(chunkStart + i, e);
                return; // the load fails as a whole; the first bad row is enough
            }
        }
        recordChunkFailure(chunkStart, batchError); // replay did not reproduce it — keep the load a failure
    }

    /** Records a statement-level error (index -1) naming the chunk range when a row cannot be attributed. */
    private void recordChunkFailure(int chunkStart, BatchUpdateException e) {
        errorCount++;
        if (errors.size() < MAX_ERRORS) {
            RowError error = SqlStateMapper.toRowError(provider, -1, e);
            error.message = error.message + " (in rows " + chunkStart + ".." + (chunkStart + batch.size() - 1) + ")";
            errors.add(error);
        }
    }

    /**
     * Fallback for drivers without savepoints: pinpoint the failing row from {@code getUpdateCounts()} when
     * exactly one EXECUTE_FAILED is reported, else fall back to a statement-level error (index -1) naming the
     * chunk range — the driver could not attribute the failure to a specific row.
     */
    private void recordFailingRowFromCounts(int chunkStart, BatchUpdateException e) {
        int[] counts = e.getUpdateCounts();
        int failPos = -1;
        boolean determinable = true;
        for (int i = 0; i < counts.length && determinable; i++)
            if (counts[i] == Statement.EXECUTE_FAILED) {
                if (failPos == -1)
                    failPos = i;
                else
                    determinable = false;
            }
        if (determinable && failPos >= 0 && failPos < batch.size())
            addError(chunkStart + failPos, e);
        else
            recordChunkFailure(chunkStart, e);
    }

    private void addError(int rowIndex, SQLException e) {
        errorCount++;
        if (errors.size() < MAX_ERRORS)
            errors.add(SqlStateMapper.toRowError(provider, rowIndex, e));
    }

    private void addMissing(int rowIndex) {
        errorCount++;
        if (errors.size() < MAX_ERRORS) {
            RowError error = new RowError();
            error.index = rowIndex;
            error.code = "missing";
            error.message = "No row matched the key columns";
            errors.add(error);
        }
    }

    private MutationResult buildResult() {
        MutationResult result = new MutationResult();
        boolean rolledBack = allOrNothing && errorCount > 0;
        result.affectedRows = rolledBack ? 0 : affected;
        result.errorCount = errorCount;
        result.errors = errors.isEmpty() ? null : errors;
        if (mode.equals("insert")) {
            result.inserted = rolledBack ? 0 : inserted;
            result.skipped = rolledBack ? 0 : skipped;
        }
        else if (mode.equals("update"))
            result.updated = rolledBack ? 0 : updated;
        // upsert reports affectedRows only — the inserted/updated split is not distinguishable in phase A
        return result;
    }

    private void releaseQuietly(Savepoint savepoint) {
        try {
            connection.releaseSavepoint(savepoint);
        } catch (SQLException e) {
            LOGGER.debug("Failed to release savepoint", e);
        }
    }

    private void closeQuietly() {
        try {
            statement.close();
        } catch (SQLException e) {
            LOGGER.warn("Failed to close bulk insert statement", e);
        }
    }
}
