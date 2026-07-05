package grok_connect.providers;

import java.nio.charset.StandardCharsets;
import java.sql.BatchUpdateException;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.Provider;
import grok_connect.table_mutation.BatchInsertBulkLoader;
import grok_connect.table_mutation.InsertRows;
import grok_connect.table_mutation.MutationManager;
import grok_connect.table_mutation.MutationResult;
import grok_connect.table_mutation.MutationValidationException;
import grok_connect.table_mutation.UpsertRows;
import grok_connect.utils.ProviderManager;
import org.mockito.Mockito;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Integration tests for GROK-20337: bulk CSV transport driven through {@link MutationManager} against
 * a containerized Postgres, for both the COPY fast path and the default chunked-executeBatch loader.
 * Proves 100k rows land identically (incl. nulls, quotes, unicode, embedded newlines), that abort and
 * a malformed stream leave zero rows, and logs COPY-vs-batch timing as flagship evidence.
 */
class PostgresBulkMutationTest extends ContainerizedProviderBaseTest {
    private static final Logger LOGGER = LoggerFactory.getLogger(PostgresBulkMutationTest.class);
    private static final int ROWS = 100_000;
    private static final int CHUNK_BYTES = 64 * 1024;
    private static final List<String> COLUMNS = Arrays.asList("id", "big", "val", "active", "note");
    private static final List<String> TYPES = Arrays.asList("int", "bigint", "double", "bool", "string");

    // Special notes at low ids exercise every CSV escaping case; nulls, unicode, quotes, comma, newline.
    private static final String UNICODE = "café ☕ → ünîcödé";
    private static final String QUOTED = "he said \"hi\"";
    private static final String NEWLINE = "line1\nline2";
    private static final String COMMA = "a,b,c";

    protected PostgresBulkMutationTest() {
        super(Provider.POSTGRESQL);
    }

    @BeforeAll
    public void initProviderManager() {
        // MutationManager resolves its provider via the GrokConnect static registry.
        GrokConnect.providerManager = new ProviderManager();
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

    private Object[] rowById(String table, int id) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             PreparedStatement statement = c.prepareStatement("SELECT big, val, active, note FROM " + table + " WHERE id = " + id);
             ResultSet rs = statement.executeQuery()) {
            Assertions.assertTrue(rs.next());
            return new Object[] {rs.getLong(1), rs.getDouble(2), rs.getBoolean(3), rs.getString(4)};
        }
    }

    @AfterAll
    public void dropScratchTables() throws SQLException {
        for (String table : new String[] {"mut_bulk_copy", "mut_bulk_batch", "mut_bulk_abort", "mut_bulk_bad",
                "mut_aon", "mut_partial", "mut_bulk_ups", "mut_upd_key", "mut_dup", "mut_blank_copy", "mut_blank_batch",
                "mut_nondeterm", "mut_ups_atomic"})
            execDirect("DROP TABLE IF EXISTS " + table);
    }

    private void createTable(String table) throws SQLException {
        execDirect("DROP TABLE IF EXISTS " + table);
        execDirect("CREATE TABLE " + table + " (id int PRIMARY KEY, big bigint, val double precision, active boolean, note text)");
    }

    private String note(int id) {
        switch (id) {
            case 0: return UNICODE;
            case 1: return QUOTED;
            case 2: return NEWLINE;
            case 3: return COMMA;
            case 4: return null; // SQL NULL
            default: return "row" + id;
        }
    }

    private String csvField(String value) {
        if (value == null)
            return ""; // empty unquoted field = COPY csv NULL
        if (value.indexOf(',') >= 0 || value.indexOf('"') >= 0 || value.indexOf('\n') >= 0 || value.indexOf('\r') >= 0)
            return "\"" + value.replace("\"", "\"\"") + "\"";
        return value;
    }

    /** Builds the full CSV payload for {@link #ROWS} rows (UTF-8, no header). */
    private byte[] buildCsv() {
        StringBuilder sb = new StringBuilder(ROWS * 24);
        for (int id = 0; id < ROWS; id++) {
            sb.append(id).append(',')
              .append(9007199254740995L + id).append(',') // beyond 2^53
              .append(id * 0.25).append(',')
              .append(id % 2 == 0).append(',')
              .append(csvField(note(id))).append('\n');
        }
        return sb.toString().getBytes(StandardCharsets.UTF_8);
    }

    private List<byte[]> chunk(byte[] data) {
        List<byte[]> chunks = new ArrayList<>();
        for (int off = 0; off < data.length; off += CHUNK_BYTES) {
            int len = Math.min(CHUNK_BYTES, data.length - off);
            byte[] c = new byte[len];
            System.arraycopy(data, off, c, 0, len);
            chunks.add(c);
        }
        return chunks;
    }

    private InsertRows bulkInsert(String table) {
        InsertRows m = new InsertRows();
        m.tableName = table;
        m.columns = COLUMNS;
        m.columnTypes = TYPES;
        m.bulk = true;
        m.rows = null; // streamed
        m.connection = connection;
        return m;
    }

    private String header(InsertRows m) {
        FuncCall call = new FuncCall();
        call.func = m;
        call.options = new HashMap<>();
        return GrokConnect.gson.toJson(call);
    }

    private void assertSpecialsAndCount(String table) throws SQLException {
        Assertions.assertEquals(ROWS, countDirect(table));
        Assertions.assertEquals(UNICODE, rowById(table, 0)[3]);
        Assertions.assertEquals(QUOTED, rowById(table, 1)[3]);
        Assertions.assertEquals(NEWLINE, rowById(table, 2)[3]);
        Assertions.assertEquals(COMMA, rowById(table, 3)[3]);
        Assertions.assertNull(rowById(table, 4)[3]); // null note
        Object[] last = rowById(table, ROWS - 1);
        Assertions.assertEquals(9007199254740995L + (ROWS - 1), last[0]); // bigint beyond 2^53, exact
        Assertions.assertEquals((ROWS - 1) * 0.25, last[1]);
        Assertions.assertEquals((ROWS - 1) % 2 == 0, last[2]);
    }

    @DisplayName("COPY fast path: 100k rows via MutationManager land with nulls/quotes/unicode/newlines intact")
    @Test
    public void copyLoader_100k() throws Exception {
        createTable("mut_bulk_copy");
        byte[] csv = buildCsv();
        MutationManager manager = new MutationManager(header(bulkInsert("mut_bulk_copy")));
        long t0 = System.nanoTime();
        manager.start();
        for (byte[] c : chunk(csv))
            manager.feed(c);
        MutationResult result = manager.finish();
        long ms = (System.nanoTime() - t0) / 1_000_000;
        LOGGER.info("BULK COPY: {} rows in {} ms", ROWS, ms);
        Assertions.assertEquals(ROWS, result.affectedRows);
        assertSpecialsAndCount("mut_bulk_copy");
    }

    @DisplayName("Default batch loader: same 100k CSV lands identically to COPY")
    @Test
    public void batchLoader_100k() throws Exception {
        createTable("mut_bulk_batch");
        byte[] csv = buildCsv();
        InsertRows m = bulkInsert("mut_bulk_batch");
        // Force the default loader (Postgres would otherwise pick COPY) by constructing it directly,
        // driven through the same connection/commit discipline MutationManager applies.
        long t0 = System.nanoTime();
        Connection conn = provider.getConnection(connection);
        try {
            provider.configureAutoCommit(conn);
            BatchInsertBulkLoader loader = new BatchInsertBulkLoader(provider, conn, m);
            MutationResult result;
            try {
                for (byte[] c : chunk(csv))
                    loader.feed(c);
                result = loader.finish();
                conn.commit();
            } catch (Exception e) {
                loader.abort();
                provider.rollbackQuietly(conn);
                throw e;
            }
            long ms = (System.nanoTime() - t0) / 1_000_000;
            LOGGER.info("BULK BATCH: {} rows in {} ms (affected {})", ROWS, ms, result.affectedRows);
        } finally {
            conn.close();
        }
        assertSpecialsAndCount("mut_bulk_batch");
    }

    @DisplayName("Abort after 50 chunks rolls back to zero rows")
    @Test
    public void abort_leavesZeroRows() throws Exception {
        createTable("mut_bulk_abort");
        List<byte[]> chunks = chunk(buildCsv());
        MutationManager manager = new MutationManager(header(bulkInsert("mut_bulk_abort")));
        manager.start();
        for (int i = 0; i < 50 && i < chunks.size(); i++)
            manager.feed(chunks.get(i));
        manager.abort();
        Assertions.assertEquals(0, countDirect("mut_bulk_abort"));
    }

    @DisplayName("Malformed CSV mid-stream surfaces an error and rolls back to zero rows")
    @Test
    public void malformed_leavesZeroRows() throws Exception {
        createTable("mut_bulk_bad");
        MutationManager manager = new MutationManager(header(bulkInsert("mut_bulk_bad")));
        manager.start();
        // valid opening rows, then a row with the wrong number of fields
        Assertions.assertThrows(Exception.class, () -> {
            manager.feed("1,100,0.5,true,ok\n2,200,1.5,false,fine\n".getBytes(StandardCharsets.UTF_8));
            manager.feed("3,only,two\n".getBytes(StandardCharsets.UTF_8)); // malformed
            manager.finish();
        });
        manager.abort();
        Assertions.assertEquals(0, countDirect("mut_bulk_bad"));
    }

    // ---- WO-6: batch modes, allOrNothing, savepoint partial errors, SqlStateMapper, carried findings ----

    private InsertRows insert(String table, List<String> cols, List<String> types) {
        InsertRows m = new InsertRows();
        m.tableName = table;
        m.columns = cols;
        m.columnTypes = types;
        m.bulk = true;
        m.rows = null; // streamed
        m.connection = connection;
        return m;
    }

    /** Builds a UTF-8, header-less CSV payload; null cells become empty unquoted fields (COPY csv NULL). */
    private byte[] csv(String[]... rows) {
        StringBuilder sb = new StringBuilder();
        for (String[] row : rows) {
            for (int i = 0; i < row.length; i++) {
                if (i > 0)
                    sb.append(',');
                sb.append(csvField(row[i]));
            }
            sb.append('\n');
        }
        return sb.toString().getBytes(StandardCharsets.UTF_8);
    }

    /** Drives a streamed mutation end-to-end through MutationManager (the real WS path minus the socket). */
    private MutationResult runManager(InsertRows m, byte[] csv) throws Exception {
        MutationManager manager = new MutationManager(header(m));
        manager.start();
        for (byte[] c : chunk(csv))
            manager.feed(c);
        return manager.finish();
    }

    /** Forces the default batch loader (Postgres would otherwise pick COPY), mirroring the manager's commit/rollback. */
    private MutationResult runBatchLoader(InsertRows m, byte[] csv) throws Exception {
        Connection conn = provider.getConnection(connection);
        try {
            provider.configureAutoCommit(conn);
            BatchInsertBulkLoader loader = new BatchInsertBulkLoader(provider, conn, m);
            MutationResult result;
            try {
                for (byte[] c : chunk(csv))
                    loader.feed(c);
                result = loader.finish();
            } catch (Exception e) {
                loader.abort();
                provider.rollbackQuietly(conn);
                throw e;
            }
            if (m.allOrNothing && result.errorCount != null && result.errorCount > 0)
                provider.rollbackQuietly(conn);
            else
                conn.commit();
            return result;
        } finally {
            conn.close();
        }
    }

    @DisplayName("(a) allOrNothing: one NOT NULL violation among 1000 rows -> 0 written, exact index + code + column")
    @Test
    public void allOrNothing_reportsFailingRowAndRollsBack() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_aon");
        execDirect("CREATE TABLE mut_aon (id int PRIMARY KEY, note text NOT NULL)");
        String[][] rows = new String[1000][];
        for (int i = 0; i < 1000; i++)
            rows[i] = new String[] {String.valueOf(i), i == 500 ? null : "n" + i};
        // insert+allOrNothing routes to COPY, which cannot report the row index — drive the batch loader,
        // the component that owns precise all-or-nothing error reporting (batchLoader_100k precedent).
        MutationResult result = runBatchLoader(insert("mut_aon", Arrays.asList("id", "note"),
                Arrays.asList("int", "string")), csv(rows));
        Assertions.assertEquals(0, result.affectedRows);
        Assertions.assertEquals(0, countDirect("mut_aon")); // whole load rolled back
        Assertions.assertEquals(Integer.valueOf(1), result.errorCount);
        Assertions.assertEquals(1, result.errors.size());
        Assertions.assertEquals(500, result.errors.get(0).index);
        Assertions.assertEquals("notnull", result.errors.get(0).code);
        Assertions.assertEquals("note", result.errors.get(0).column);
    }

    @DisplayName("(b) partial: 999 land, the 1 bad row is reported, counters exact")
    @Test
    public void partial_keepsSurvivorsAndReportsBadRow() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_partial");
        execDirect("CREATE TABLE mut_partial (id int PRIMARY KEY, note text NOT NULL)");
        String[][] rows = new String[1000][];
        for (int i = 0; i < 1000; i++)
            rows[i] = new String[] {String.valueOf(i), i == 500 ? null : "n" + i};
        InsertRows m = insert("mut_partial", Arrays.asList("id", "note"), Arrays.asList("int", "string"));
        m.allOrNothing = false;
        MutationResult result = runManager(m, csv(rows));
        Assertions.assertEquals(999, result.affectedRows);
        Assertions.assertEquals(999, countDirect("mut_partial"));
        Assertions.assertEquals(Integer.valueOf(999), result.inserted);
        Assertions.assertEquals(Integer.valueOf(1), result.errorCount);
        Assertions.assertEquals(500, result.errors.get(0).index);
        Assertions.assertEquals("notnull", result.errors.get(0).code);
        Assertions.assertEquals(0, ((Number) queryCount("SELECT count(*) FROM mut_partial WHERE id = 500")).intValue());
    }

    @DisplayName("(c) upsert via bulk: 3 seeded, 5-row payload (2 update + 3 insert) -> final state correct, affectedRows=5")
    @Test
    public void upsert_viaBulk() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_bulk_ups");
        execDirect("CREATE TABLE mut_bulk_ups (id int PRIMARY KEY, region text, amount double precision)");
        execDirect("INSERT INTO mut_bulk_ups VALUES (1, 'east', 10), (2, 'west', 20), (3, 'north', 30)");
        UpsertRows m = new UpsertRows();
        m.tableName = "mut_bulk_ups";
        m.columns = Arrays.asList("id", "region", "amount");
        m.columnTypes = Arrays.asList("int", "string", "double");
        m.matchKeys = Arrays.asList("id");
        m.bulk = true;
        m.rows = null;
        m.connection = connection;
        MutationResult result = runManager(m, csv(
                new String[] {"1", "east2", "111"},   // update
                new String[] {"2", "west2", "222"},   // update
                new String[] {"4", "south", "40"},    // insert
                new String[] {"5", "central", "50"},  // insert
                new String[] {"6", "east3", "60"}));   // insert
        Assertions.assertEquals(5, result.affectedRows);
        Assertions.assertEquals(6, countDirect("mut_bulk_ups"));
        Assertions.assertEquals(111.0d, scalarDouble("SELECT amount FROM mut_bulk_ups WHERE id = 1"));
        Assertions.assertEquals(222.0d, scalarDouble("SELECT amount FROM mut_bulk_ups WHERE id = 2"));
        Assertions.assertEquals(30.0d, scalarDouble("SELECT amount FROM mut_bulk_ups WHERE id = 3")); // untouched
        Assertions.assertEquals(40.0d, scalarDouble("SELECT amount FROM mut_bulk_ups WHERE id = 4"));
    }

    @DisplayName("(d) update by key: matching rows land, a missing key is reported as code=missing")
    @Test
    public void updateByKey_reportsMissing() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_upd_key");
        execDirect("CREATE TABLE mut_upd_key (id int PRIMARY KEY, amount double precision)");
        execDirect("INSERT INTO mut_upd_key VALUES (1, 10), (2, 20), (3, 30)");
        InsertRows m = insert("mut_upd_key", Arrays.asList("id", "amount"), Arrays.asList("int", "double"));
        m.mode = "update";
        m.keyColumns = Arrays.asList("id");
        m.allOrNothing = false;
        MutationResult result = runManager(m, csv(
                new String[] {"1", "100"},
                new String[] {"2", "200"},
                new String[] {"99", "999"})); // no such row
        Assertions.assertEquals(2, result.affectedRows);
        Assertions.assertEquals(Integer.valueOf(2), result.updated);
        Assertions.assertEquals(Integer.valueOf(1), result.errorCount);
        Assertions.assertEquals(2, result.errors.get(0).index);
        Assertions.assertEquals("missing", result.errors.get(0).code);
        Assertions.assertEquals(100.0d, scalarDouble("SELECT amount FROM mut_upd_key WHERE id = 1"));
        Assertions.assertEquals(200.0d, scalarDouble("SELECT amount FROM mut_upd_key WHERE id = 2"));
        Assertions.assertEquals(30.0d, scalarDouble("SELECT amount FROM mut_upd_key WHERE id = 3")); // untouched
        Assertions.assertEquals(3, countDirect("mut_upd_key")); // missing key not inserted
    }

    @DisplayName("(e) errorOnDuplicate=false skips the duplicate; =true reports it as a unique error")
    @Test
    public void errorOnDuplicate_bothSettings() throws Exception {
        // false: ON CONFLICT DO NOTHING -> the duplicate is skipped, the rest land
        execDirect("DROP TABLE IF EXISTS mut_dup");
        execDirect("CREATE TABLE mut_dup (id int PRIMARY KEY, note text)");
        execDirect("INSERT INTO mut_dup VALUES (1, 'seed')");
        InsertRows skip = insert("mut_dup", Arrays.asList("id", "note"), Arrays.asList("int", "string"));
        skip.allOrNothing = false;
        skip.errorOnDuplicate = false;
        MutationResult skipResult = runManager(skip, csv(
                new String[] {"1", "a"}, new String[] {"2", "b"}, new String[] {"3", "c"}));
        Assertions.assertEquals(Integer.valueOf(2), skipResult.inserted);
        Assertions.assertEquals(Integer.valueOf(1), skipResult.skipped);
        Assertions.assertEquals(Integer.valueOf(0), skipResult.errorCount);
        Assertions.assertEquals(3, countDirect("mut_dup"));
        Assertions.assertEquals("seed", scalarString("SELECT note FROM mut_dup WHERE id = 1")); // not overwritten

        // true: the duplicate surfaces as a unique RowError, the rest still land (partial)
        execDirect("DROP TABLE IF EXISTS mut_dup");
        execDirect("CREATE TABLE mut_dup (id int PRIMARY KEY, note text)");
        execDirect("INSERT INTO mut_dup VALUES (1, 'seed')");
        InsertRows err = insert("mut_dup", Arrays.asList("id", "note"), Arrays.asList("int", "string"));
        err.allOrNothing = false;
        err.errorOnDuplicate = true;
        MutationResult errResult = runManager(err, csv(
                new String[] {"1", "a"}, new String[] {"2", "b"}, new String[] {"3", "c"}));
        Assertions.assertEquals(Integer.valueOf(2), errResult.inserted);
        Assertions.assertEquals(Integer.valueOf(1), errResult.errorCount);
        Assertions.assertEquals(0, errResult.errors.get(0).index);
        Assertions.assertEquals("unique", errResult.errors.get(0).code);
        Assertions.assertEquals(3, countDirect("mut_dup"));
    }

    @DisplayName("(f) a provider without savepoint support rejects partial mode with a capability error")
    @Test
    public void partialMode_savepointless_capabilityError() throws Exception {
        DatabaseMetaData meta = Mockito.mock(DatabaseMetaData.class);
        Mockito.when(meta.supportsSavepoints()).thenReturn(false);
        Connection conn = Mockito.mock(Connection.class);
        Mockito.when(conn.getMetaData()).thenReturn(meta);
        InsertRows m = insert("mut_nosp", Arrays.asList("id"), Arrays.asList("int"));
        m.allOrNothing = false;
        UnsupportedOperationException ex = Assertions.assertThrows(UnsupportedOperationException.class,
                () -> new BatchInsertBulkLoader(provider, conn, m));
        Assertions.assertTrue(ex.getMessage().contains("Partial mode"));
    }

    @DisplayName("Carried finding 1: a blank line / null row lands identically through COPY and the batch loader")
    @Test
    public void blankLine_copyAndBatchAgree() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_blank_copy");
        execDirect("DROP TABLE IF EXISTS mut_blank_batch");
        execDirect("CREATE TABLE mut_blank_copy (v text)");
        execDirect("CREATE TABLE mut_blank_batch (v text)");
        byte[] payload = "a\n\nb\n".getBytes(StandardCharsets.UTF_8); // middle line is blank = a null row

        runManager(insert("mut_blank_copy", Arrays.asList("v"), Arrays.asList("string")), payload); // COPY
        runBatchLoader(insert("mut_blank_batch", Arrays.asList("v"), Arrays.asList("string")), payload); // batch

        Assertions.assertEquals(3, countDirect("mut_blank_copy"));
        Assertions.assertEquals(countDirect("mut_blank_copy"), countDirect("mut_blank_batch")); // no divergence
        Assertions.assertEquals(1, ((Number) queryCount("SELECT count(*) FROM mut_blank_copy WHERE v IS NULL")).intValue());
        Assertions.assertEquals(1, ((Number) queryCount("SELECT count(*) FROM mut_blank_batch WHERE v IS NULL")).intValue());
    }

    @DisplayName("Carried finding 2: an unwired bulk mode fails loud instead of silently inserting")
    @Test
    public void unknownBulkMode_failsLoud() {
        InsertRows m = insert("mut_unknown", Arrays.asList("id"), Arrays.asList("int"));
        m.mode = "frobnicate";
        Assertions.assertThrows(MutationValidationException.class, () -> {
            MutationManager manager = new MutationManager(header(m));
            manager.start();
        });
    }

    @DisplayName("Fix 1: a non-reproducing atomic batch failure still forces full rollback (no silent partial commit)")
    @Test
    public void atomicReplay_nonReproducingFailure_rollsBackFully() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_nondeterm");
        execDirect("CREATE TABLE mut_nondeterm (id int PRIMARY KEY, note text)");
        // Raw pgjdbc connection (HikariProxyConnection is final and cannot be spied); drive it like the manager.
        Connection real = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
        try {
            real.setAutoCommit(false);
            Connection spy = Mockito.spy(real);
            // The batch throws a transient 40P01 (deadlock) that does NOT recur on the row-by-row replay —
            // executeBatch is stubbed to throw, executeUpdate delegates to the real statement and succeeds.
            Mockito.doAnswer(inv -> {
                PreparedStatement realStmt = (PreparedStatement) inv.callRealMethod();
                PreparedStatement spyStmt = Mockito.spy(realStmt);
                Mockito.doThrow(new BatchUpdateException("deadlock detected", "40P01", new int[0]))
                        .when(spyStmt).executeBatch();
                return spyStmt;
            }).when(spy).prepareStatement(Mockito.anyString());

            InsertRows m = insert("mut_nondeterm", Arrays.asList("id", "note"), Arrays.asList("int", "string"));
            m.errorOnDuplicate = true; // plain insert, no ON CONFLICT — the replay would otherwise succeed silently
            BatchInsertBulkLoader loader = new BatchInsertBulkLoader(provider, spy, m);
            byte[] csv = csv(new String[] {"1", "a"}, new String[] {"2", "b"}, new String[] {"3", "c"});
            MutationResult result;
            try {
                loader.feed(csv);
                result = loader.finish();
            } catch (Exception e) {
                loader.abort();
                real.rollback();
                throw e;
            }
            // errorCount>0 must force the whole load to roll back even though no single row reproduced the error
            Assertions.assertTrue(result.errorCount != null && result.errorCount > 0);
            Assertions.assertEquals(0, result.affectedRows);
            Assertions.assertEquals("conflict", result.errors.get(0).code); // 40P01 -> conflict
            real.rollback(); // the manager would roll back; assert nothing persisted
        } finally {
            real.close();
        }
        Assertions.assertEquals(0, countDirect("mut_nondeterm"));
    }

    @DisplayName("Fix 2: allOrNothing upsert with a mid-batch violation rolls back fully through MutationManager")
    @Test
    public void upsert_allOrNothing_atomicRollbackThroughManager() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_ups_atomic");
        execDirect("CREATE TABLE mut_ups_atomic (id int PRIMARY KEY, note text NOT NULL)");
        UpsertRows m = new UpsertRows(); // allOrNothing defaults to true
        m.tableName = "mut_ups_atomic";
        m.columns = Arrays.asList("id", "note");
        m.columnTypes = Arrays.asList("int", "string");
        m.matchKeys = Arrays.asList("id");
        m.bulk = true;
        m.rows = null;
        m.connection = connection;
        // upsert routes to the batch loader (not COPY), so this exercises the manager's atomic rollback end-to-end
        MutationResult result = runManager(m, csv(
                new String[] {"1", "a"},
                new String[] {"2", null}, // NOT NULL violation, mid-batch
                new String[] {"3", "c"}));
        Assertions.assertEquals(0, result.affectedRows);
        Assertions.assertEquals(0, countDirect("mut_ups_atomic")); // whole load rolled back
        Assertions.assertEquals(Integer.valueOf(1), result.errorCount);
        Assertions.assertEquals(1, result.errors.get(0).index); // exact failing row
        Assertions.assertEquals("notnull", result.errors.get(0).code);
    }

    private Object queryCount(String sql) throws SQLException {
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             ResultSet rs = c.createStatement().executeQuery(sql)) {
            rs.next();
            return rs.getObject(1);
        }
    }

    private double scalarDouble(String sql) throws SQLException {
        return ((Number) queryCount(sql)).doubleValue();
    }

    private String scalarString(String sql) throws SQLException {
        return (String) queryCount(sql);
    }
}
