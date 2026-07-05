package grok_connect.providers;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
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
import java.util.Collections;
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
import serialization.BigIntColumn;
import serialization.BoolColumn;
import serialization.Column;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;

/**
 * Integration tests for GROK-20347 (java-d42-reader WO-4): the bulk mutation transport re-driven with
 * typed d42 frames, through {@link MutationManager} against a containerized Postgres, for both the COPY
 * fast path ({@code CopyCsvFormatter}) and the default typed-batch loader. Re-drives the WO-6 batch
 * matrix (allOrNothing / partial / upsert / update / errorOnDuplicate / atomic replay), proves the two
 * loaders store identical values (loader-parity pin), replays committed Dart d42 fixtures (exotic
 * encoders) through the real transport, and keeps the legacy CSV path green through the batch loader.
 */
class PostgresBulkMutationTest extends ContainerizedProviderBaseTest {
    private static final Logger LOGGER = LoggerFactory.getLogger(PostgresBulkMutationTest.class);
    private static final int ROWS = 100_000;
    private static final int CHUNK_ROWS = 20_000;
    private static final int CHUNK_BYTES = 64 * 1024;
    private static final List<String> COLUMNS = Arrays.asList("id", "big", "val", "active", "note");
    private static final List<String> TYPES = Arrays.asList("int", "bigint", "double", "bool", "string");

    private static final String UNICODE = "café ☕ → ünîcödé";
    private static final String QUOTED = "he said \"hi\"";
    private static final String NEWLINE = "line1\nline2";
    private static final String COMMA = "a,b,c";

    protected PostgresBulkMutationTest() {
        super(Provider.POSTGRESQL);
    }

    @BeforeAll
    public void initProviderManager() {
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
                "mut_aon", "mut_partial", "mut_bulk_ups", "mut_upd_key", "mut_dup", "mut_null_copy", "mut_null_batch",
                "mut_nondeterm", "mut_ups_atomic", "mut_parity_copy", "mut_parity_batch", "mut_csv_batch", "mut_mismatch",
                "mut_fx_forced_string_prefixes", "mut_fx_forced_string_squash",
                "mut_fx_natural_string_categories", "mut_fx_natural_float64"})
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

    // ---- d42 frame construction (the transport now carries typed DataFrames, not CSV) ----

    /** Builds a typed d42 column from string cells (null cell -> None), by dg type. */
    private Column<?> column(String name, String type, String[] rows, int col) {
        int n = rows == null ? 0 : rows.length;
        switch (type) {
            case "int": {
                Integer[] v = new Integer[n];
                for (int r = 0; r < n; r++) v[r] = cell(rows, r) == null ? null : Integer.parseInt(cell(rows, r));
                return new IntColumn(name, v);
            }
            case "bigint": {
                String[] v = new String[n];
                for (int r = 0; r < n; r++) v[r] = cell(rows, r);
                return new BigIntColumn(name, v);
            }
            case "double": {
                Double[] v = new Double[n];
                for (int r = 0; r < n; r++) v[r] = cell(rows, r) == null ? null : Double.parseDouble(cell(rows, r));
                return FloatColumn.double64(name, v);
            }
            case "bool": {
                Boolean[] v = new Boolean[n];
                for (int r = 0; r < n; r++) v[r] = cell(rows, r) == null ? null : Boolean.parseBoolean(cell(rows, r));
                return new BoolColumn(name, v);
            }
            case "datetime": {
                Double[] v = new Double[n];
                for (int r = 0; r < n; r++) v[r] = cell(rows, r) == null ? null : Double.parseDouble(cell(rows, r));
                return new DateTimeColumn(name, v);
            }
            default: {
                String[] v = new String[n];
                for (int r = 0; r < n; r++) v[r] = cell(rows, r);
                return new StringColumn(name, v);
            }
        }
    }

    private String cell(String[] flat, int r) { return flat[r]; }

    /** Serializes a typed frame (one column per name/type, values from row-major {@code rows}) to d42. */
    private byte[] d42(List<String> names, List<String> types, String[]... rows) {
        Column<?>[] cols = new Column[names.size()];
        for (int c = 0; c < names.size(); c++) {
            String[] colValues = new String[rows.length];
            for (int r = 0; r < rows.length; r++)
                colValues[r] = rows[r][c];
            cols[c] = column(names.get(c), types.get(c), colValues, c);
        }
        return DataFrame.fromColumns(cols).toByteArray();
    }

    /** Builds the {@link #ROWS}-row flagship frame as d42 sub-frames of {@link #CHUNK_ROWS} rows each. */
    private List<byte[]> flagshipChunks() {
        List<byte[]> chunks = new ArrayList<>();
        for (int start = 0; start < ROWS; start += CHUNK_ROWS) {
            int len = Math.min(CHUNK_ROWS, ROWS - start);
            Integer[] id = new Integer[len];
            String[] big = new String[len];
            Double[] val = new Double[len];
            Boolean[] active = new Boolean[len];
            String[] note = new String[len];
            for (int i = 0; i < len; i++) {
                int gid = start + i;
                id[i] = gid;
                big[i] = Long.toString(9007199254740995L + gid); // beyond 2^53
                val[i] = gid * 0.25;
                active[i] = gid % 2 == 0;
                note[i] = note(gid);
            }
            chunks.add(DataFrame.fromColumns(
                    new IntColumn("id", id),
                    new BigIntColumn("big", big),
                    FloatColumn.double64("val", val),
                    new BoolColumn("active", active),
                    new StringColumn("note", note)).toByteArray());
        }
        return chunks;
    }

    private InsertRows bulkInsert(String table) {
        return insert(table, COLUMNS, TYPES);
    }

    private InsertRows insert(String table, List<String> cols, List<String> types) {
        InsertRows m = new InsertRows();
        m.tableName = table;
        m.columns = cols;
        m.columnTypes = types;
        m.bulk = true;
        m.rows = null; // streamed
        m.payloadFormat = "d42";
        m.connection = connection;
        return m;
    }

    private String header(InsertRows m) {
        FuncCall call = new FuncCall();
        call.func = m;
        call.options = new HashMap<>();
        return GrokConnect.gson.toJson(call);
    }

    /** Drives a streamed d42 mutation end-to-end through MutationManager (the real WS path minus the socket). */
    private MutationResult runManager(InsertRows m, byte[]... chunks) throws Exception {
        MutationManager manager = new MutationManager(header(m));
        manager.start();
        for (byte[] c : chunks)
            manager.feed(c);
        return manager.finish();
    }

    /** Forces the default typed-batch loader (Postgres would otherwise pick COPY), mirroring the manager's commit/rollback. */
    private MutationResult runBatchLoader(InsertRows m, byte[]... chunks) throws Exception {
        Connection conn = provider.getConnection(connection);
        try {
            provider.configureAutoCommit(conn);
            BatchInsertBulkLoader loader = new BatchInsertBulkLoader(provider, conn, m);
            MutationResult result;
            try {
                for (byte[] c : chunks)
                    loader.feed(DataFrame.fromByteArray(c));
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

    @DisplayName("COPY fast path: 100k-row d42 via MutationManager land with nulls/quotes/unicode/newlines intact")
    @Test
    public void copyLoader_100k() throws Exception {
        createTable("mut_bulk_copy");
        List<byte[]> chunks = flagshipChunks();
        MutationManager manager = new MutationManager(header(bulkInsert("mut_bulk_copy")));
        long t0 = System.nanoTime();
        manager.start();
        for (byte[] c : chunks)
            manager.feed(c);
        MutationResult result = manager.finish();
        long ms = (System.nanoTime() - t0) / 1_000_000;
        LOGGER.info("BULK d42 COPY: {} rows in {} ms", ROWS, ms);
        Assertions.assertEquals(ROWS, result.affectedRows);
        assertSpecialsAndCount("mut_bulk_copy");
    }

    @DisplayName("Default typed-batch loader: same 100k d42 lands identically to COPY")
    @Test
    public void batchLoader_100k() throws Exception {
        createTable("mut_bulk_batch");
        InsertRows m = bulkInsert("mut_bulk_batch");
        long t0 = System.nanoTime();
        Connection conn = provider.getConnection(connection);
        try {
            provider.configureAutoCommit(conn);
            BatchInsertBulkLoader loader = new BatchInsertBulkLoader(provider, conn, m);
            MutationResult result;
            try {
                for (byte[] c : flagshipChunks())
                    loader.feed(DataFrame.fromByteArray(c));
                result = loader.finish();
                conn.commit();
            } catch (Exception e) {
                loader.abort();
                provider.rollbackQuietly(conn);
                throw e;
            }
            long ms = (System.nanoTime() - t0) / 1_000_000;
            LOGGER.info("BULK d42 BATCH: {} rows in {} ms (affected {})", ROWS, ms, result.affectedRows);
        } finally {
            conn.close();
        }
        assertSpecialsAndCount("mut_bulk_batch");
    }

    @DisplayName("Legacy CSV path (dual-format): 100k CSV via the batch loader feedCsv still lands (+ timing)")
    @Test
    public void csvBatch_100k_dualFormat() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_csv_batch");
        execDirect("CREATE TABLE mut_csv_batch (id int PRIMARY KEY, big bigint, val double precision, active boolean, note text)");
        InsertRows m = bulkInsert("mut_csv_batch");
        m.payloadFormat = "csv"; // exercise the surviving CSV branch through the batch loader
        byte[] csv = buildCsv();
        long t0 = System.nanoTime();
        Connection conn = provider.getConnection(connection);
        try {
            provider.configureAutoCommit(conn);
            BatchInsertBulkLoader loader = new BatchInsertBulkLoader(provider, conn, m);
            for (byte[] c : chunkBytes(csv))
                loader.feedCsv(c);
            MutationResult result = loader.finish();
            conn.commit();
            long ms = (System.nanoTime() - t0) / 1_000_000;
            LOGGER.info("BULK CSV BATCH (legacy): {} rows in {} ms (affected {})", ROWS, ms, result.affectedRows);
        } finally {
            conn.close();
        }
        assertSpecialsAndCount("mut_csv_batch");
    }

    @DisplayName("Abort after a chunk rolls back to zero rows")
    @Test
    public void abort_leavesZeroRows() throws Exception {
        createTable("mut_bulk_abort");
        List<byte[]> chunks = flagshipChunks();
        MutationManager manager = new MutationManager(header(bulkInsert("mut_bulk_abort")));
        manager.start();
        manager.feed(chunks.get(0));
        manager.abort();
        Assertions.assertEquals(0, countDirect("mut_bulk_abort"));
    }

    @DisplayName("Schema-mismatch d42 frame surfaces an error and rolls back to zero rows")
    @Test
    public void schemaMismatch_leavesZeroRows() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_mismatch");
        execDirect("CREATE TABLE mut_mismatch (id int PRIMARY KEY, note text)");
        MutationManager manager = new MutationManager(header(insert("mut_mismatch",
                Arrays.asList("id", "note"), Arrays.asList("int", "string"))));
        manager.start();
        // frame's 2nd column is named 'wrong' — validateChunkSchema must reject it before any insert
        byte[] bad = d42(Arrays.asList("id", "wrong"), Arrays.asList("int", "string"),
                new String[] {"1", "a"}, new String[] {"2", "b"});
        Assertions.assertThrows(MutationValidationException.class, () -> manager.feed(bad));
        manager.abort();
        Assertions.assertEquals(0, countDirect("mut_mismatch"));
    }

    @DisplayName("Loader-parity pin: identical d42 chunk through COPY and batch stores identical values (us timestamps, doubles, nulls)")
    @Test
    public void loaderParity_copyVsBatch_identical() throws Exception {
        for (String t : new String[] {"mut_parity_copy", "mut_parity_batch"}) {
            execDirect("DROP TABLE IF EXISTS " + t);
            execDirect("CREATE TABLE " + t + " (id int PRIMARY KEY, note text, amount double precision, active boolean, ts timestamp)");
        }
        List<String> cols = Arrays.asList("id", "note", "amount", "active", "ts");
        List<String> types = Arrays.asList("int", "string", "double", "bool", "datetime");
        byte[] chunk = d42(cols, types,
                new String[] {"1", "east",  "1.0000000000000002", "true",  "1751722496789012"},
                new String[] {"2", null,     "-0.0",               "false", "1751722496000000"},
                new String[] {"3", "a,b\"c", "1.0e308",            "true",  "1000000000123456"},
                new String[] {"4", "z",      null,                 "false", null});
        // COPY (insert + allOrNothing routes to CopyCsvFormatter); batch (forced)
        runManager(insert("mut_parity_copy", cols, types), chunk);
        runBatchLoader(insert("mut_parity_batch", cols, types), chunk);
        Assertions.assertEquals(4, countDirect("mut_parity_copy"));
        Assertions.assertEquals(4, countDirect("mut_parity_batch"));
        // IS DISTINCT FROM covers nulls; exact equality covers us-precision timestamps and doubles
        int diffs = ((Number) queryCount(
                "SELECT count(*) FROM mut_parity_copy c FULL JOIN mut_parity_batch b USING (id) WHERE "
                + "c.note IS DISTINCT FROM b.note OR c.amount IS DISTINCT FROM b.amount OR "
                + "c.active IS DISTINCT FROM b.active OR c.ts IS DISTINCT FROM b.ts")).intValue();
        Assertions.assertEquals(0, diffs, "COPY and batch stored different values for the same chunk");
    }

    // ---- WO-6 batch matrix, re-driven with d42 frames (semantics unchanged) ----

    @DisplayName("(a) allOrNothing: one NOT NULL violation among 1000 rows -> 0 written, exact index + code + column")
    @Test
    public void allOrNothing_reportsFailingRowAndRollsBack() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_aon");
        execDirect("CREATE TABLE mut_aon (id int PRIMARY KEY, note text NOT NULL)");
        String[][] rows = new String[1000][];
        for (int i = 0; i < 1000; i++)
            rows[i] = new String[] {String.valueOf(i), i == 500 ? null : "n" + i};
        // insert+allOrNothing routes to COPY, which cannot report the row index — drive the batch loader.
        MutationResult result = runBatchLoader(insert("mut_aon", Arrays.asList("id", "note"),
                Arrays.asList("int", "string")), d42(Arrays.asList("id", "note"), Arrays.asList("int", "string"), rows));
        Assertions.assertEquals(0, result.affectedRows);
        Assertions.assertEquals(0, countDirect("mut_aon"));
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
        MutationResult result = runManager(m, d42(Arrays.asList("id", "note"), Arrays.asList("int", "string"), rows));
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
        m.payloadFormat = "d42";
        m.connection = connection;
        MutationResult result = runManager(m, d42(m.columns, m.columnTypes,
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
        MutationResult result = runManager(m, d42(m.columns, m.columnTypes,
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
        Assertions.assertEquals(3, countDirect("mut_upd_key"));
    }

    @DisplayName("(e) errorOnDuplicate=false skips the duplicate; =true reports it as a unique error")
    @Test
    public void errorOnDuplicate_bothSettings() throws Exception {
        List<String> cols = Arrays.asList("id", "note");
        List<String> types = Arrays.asList("int", "string");
        execDirect("DROP TABLE IF EXISTS mut_dup");
        execDirect("CREATE TABLE mut_dup (id int PRIMARY KEY, note text)");
        execDirect("INSERT INTO mut_dup VALUES (1, 'seed')");
        InsertRows skip = insert("mut_dup", cols, types);
        skip.allOrNothing = false;
        skip.errorOnDuplicate = false;
        MutationResult skipResult = runManager(skip, d42(cols, types,
                new String[] {"1", "a"}, new String[] {"2", "b"}, new String[] {"3", "c"}));
        Assertions.assertEquals(Integer.valueOf(2), skipResult.inserted);
        Assertions.assertEquals(Integer.valueOf(1), skipResult.skipped);
        Assertions.assertEquals(Integer.valueOf(0), skipResult.errorCount);
        Assertions.assertEquals(3, countDirect("mut_dup"));
        Assertions.assertEquals("seed", scalarString("SELECT note FROM mut_dup WHERE id = 1"));

        execDirect("DROP TABLE IF EXISTS mut_dup");
        execDirect("CREATE TABLE mut_dup (id int PRIMARY KEY, note text)");
        execDirect("INSERT INTO mut_dup VALUES (1, 'seed')");
        InsertRows err = insert("mut_dup", cols, types);
        err.allOrNothing = false;
        err.errorOnDuplicate = true;
        MutationResult errResult = runManager(err, d42(cols, types,
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

    @DisplayName("A none string cell lands as SQL NULL identically through COPY and the batch loader")
    @Test
    public void nullCell_copyAndBatchAgree() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_null_copy");
        execDirect("DROP TABLE IF EXISTS mut_null_batch");
        execDirect("CREATE TABLE mut_null_copy (v text)");
        execDirect("CREATE TABLE mut_null_batch (v text)");
        byte[] chunk = d42(Arrays.asList("v"), Arrays.asList("string"),
                new String[] {"a"}, new String[] {(String) null}, new String[] {"b"});

        runManager(insert("mut_null_copy", Arrays.asList("v"), Arrays.asList("string")), chunk); // COPY
        runBatchLoader(insert("mut_null_batch", Arrays.asList("v"), Arrays.asList("string")), chunk); // batch

        Assertions.assertEquals(3, countDirect("mut_null_copy"));
        Assertions.assertEquals(countDirect("mut_null_copy"), countDirect("mut_null_batch"));
        Assertions.assertEquals(1, ((Number) queryCount("SELECT count(*) FROM mut_null_copy WHERE v IS NULL")).intValue());
        Assertions.assertEquals(1, ((Number) queryCount("SELECT count(*) FROM mut_null_batch WHERE v IS NULL")).intValue());
    }

    @DisplayName("An unwired bulk mode fails loud instead of silently inserting")
    @Test
    public void unknownBulkMode_failsLoud() {
        InsertRows m = insert("mut_unknown", Arrays.asList("id"), Arrays.asList("int"));
        m.mode = "frobnicate";
        Assertions.assertThrows(MutationValidationException.class, () -> {
            MutationManager manager = new MutationManager(header(m));
            manager.start();
        });
    }

    @DisplayName("Unknown payloadFormat fails at construction")
    @Test
    public void unknownPayloadFormat_failsLoud() {
        InsertRows m = insert("mut_unknown", Arrays.asList("id"), Arrays.asList("int"));
        m.payloadFormat = "avro";
        Assertions.assertThrows(MutationValidationException.class, () -> new MutationManager(header(m)));
    }

    @DisplayName("Fix 1: a non-reproducing atomic batch failure still forces full rollback (no silent partial commit)")
    @Test
    public void atomicReplay_nonReproducingFailure_rollsBackFully() throws Exception {
        execDirect("DROP TABLE IF EXISTS mut_nondeterm");
        execDirect("CREATE TABLE mut_nondeterm (id int PRIMARY KEY, note text)");
        Connection real = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
        try {
            real.setAutoCommit(false);
            Connection spy = Mockito.spy(real);
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
            byte[] chunk = d42(Arrays.asList("id", "note"), Arrays.asList("int", "string"),
                    new String[] {"1", "a"}, new String[] {"2", "b"}, new String[] {"3", "c"});
            MutationResult result;
            try {
                loader.feed(DataFrame.fromByteArray(chunk));
                result = loader.finish();
            } catch (Exception e) {
                loader.abort();
                real.rollback();
                throw e;
            }
            Assertions.assertTrue(result.errorCount != null && result.errorCount > 0);
            Assertions.assertEquals(0, result.affectedRows);
            Assertions.assertEquals("conflict", result.errors.get(0).code); // 40P01 -> conflict
            real.rollback();
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
        UpsertRows m = new UpsertRows();
        m.tableName = "mut_ups_atomic";
        m.columns = Arrays.asList("id", "note");
        m.columnTypes = Arrays.asList("int", "string");
        m.matchKeys = Arrays.asList("id");
        m.bulk = true;
        m.rows = null;
        m.payloadFormat = "d42";
        m.connection = connection;
        MutationResult result = runManager(m, d42(m.columns, m.columnTypes,
                new String[] {"1", "a"},
                new String[] {"2", null}, // NOT NULL violation, mid-batch
                new String[] {"3", "c"}));
        Assertions.assertEquals(0, result.affectedRows);
        Assertions.assertEquals(0, countDirect("mut_ups_atomic"));
        Assertions.assertEquals(Integer.valueOf(1), result.errorCount);
        Assertions.assertEquals(1, result.errors.get(0).index);
        Assertions.assertEquals("notnull", result.errors.get(0).code);
    }

    // ---- Committed Dart d42 fixtures (exotic encoders) replayed through the real transport ----

    private byte[] readD42Fixture(String name) throws IOException {
        for (String root : new String[] {
                "../serialization/src/test/resources/d42/",
                "serialization/src/test/resources/d42/",
                "connectors/serialization/src/test/resources/d42/"}) {
            File f = new File(root + name + ".d42");
            if (f.exists())
                return Files.readAllBytes(f.toPath());
        }
        throw new IOException("d42 fixture not found on any known root: " + name);
    }

    private String sqlType(String dgType) {
        switch (dgType) {
            case "int": return "int";
            case "bigint": return "bigint";
            case "double": return "double precision";
            case "bool": return "boolean";
            case "datetime": return "timestamp";
            default: return "text";
        }
    }

    /** Replays a single-column Dart fixture through the transport and asserts the DB round-trips the decoded values. */
    private void replaySingleColumn(String fixture, boolean forceBatch) throws Exception {
        byte[] bytes = readD42Fixture(fixture);
        DataFrame expected = DataFrame.fromByteArray(bytes);
        Column<?> col = expected.getColumn(0);
        String table = "mut_fx_" + fixture;
        execDirect("DROP TABLE IF EXISTS " + table);
        execDirect("CREATE TABLE " + table + " (\"" + col.getName() + "\" " + sqlType(col.getType()) + ")");
        InsertRows m = insert(table, Collections.singletonList(col.getName()), Collections.singletonList(col.getType()));
        if (forceBatch)
            runBatchLoader(m, bytes);
        else
            runManager(m, bytes);
        Assertions.assertEquals(expected.rowCount.intValue(), (int) countDirect(table),
                fixture + ": row count through transport");

        if (col.getType().equals("string")) {
            List<String> exp = new ArrayList<>();
            for (int r = 0; r < expected.rowCount; r++)
                exp.add(col.isNone(r) ? null : col.get(r).toString());
            List<String> got = queryStrings("SELECT \"" + col.getName() + "\" FROM " + table);
            exp.sort(nullsLast());
            got.sort(nullsLast());
            Assertions.assertEquals(exp, got, fixture + ": string values through transport");
        }
        else { // double
            List<Double> exp = new ArrayList<>();
            for (int r = 0; r < expected.rowCount; r++)
                exp.add(col.isNone(r) ? null : ((FloatColumn) col).getDouble(r));
            List<Double> got = queryDoubles("SELECT \"" + col.getName() + "\" FROM " + table);
            exp.sort(nullDoubles());
            got.sort(nullDoubles());
            Assertions.assertEquals(exp.size(), got.size());
            for (int i = 0; i < exp.size(); i++)
                Assertions.assertTrue(doubleEq(exp.get(i), got.get(i)),
                        fixture + ": double mismatch at " + i + " expected " + exp.get(i) + " got " + got.get(i));
        }
    }

    @DisplayName("Dart fixtures: string:prefixes + string:squash replay through the COPY transport (exotic decoders)")
    @Test
    public void dartFixtures_stringEncoders_viaCopy() throws Exception {
        replaySingleColumn("forced_string_prefixes", false);
        replaySingleColumn("forced_string_squash", false);
    }

    @DisplayName("Dart fixtures: string:categories(+int:bitIntList) + float:raw64 replay through the batch bind path")
    @Test
    public void dartFixtures_viaBatchBind() throws Exception {
        replaySingleColumn("natural_string_categories", true);
        replaySingleColumn("natural_float64", true);
    }

    // ---- CSV helpers (legacy dual-format path) ----

    private String csvField(String value) {
        if (value == null)
            return "";
        if (value.indexOf(',') >= 0 || value.indexOf('"') >= 0 || value.indexOf('\n') >= 0 || value.indexOf('\r') >= 0)
            return "\"" + value.replace("\"", "\"\"") + "\"";
        return value;
    }

    private byte[] buildCsv() {
        StringBuilder sb = new StringBuilder(ROWS * 24);
        for (int id = 0; id < ROWS; id++) {
            sb.append(id).append(',')
              .append(9007199254740995L + id).append(',')
              .append(id * 0.25).append(',')
              .append(id % 2 == 0).append(',')
              .append(csvField(note(id))).append('\n');
        }
        return sb.toString().getBytes(StandardCharsets.UTF_8);
    }

    private List<byte[]> chunkBytes(byte[] data) {
        List<byte[]> chunks = new ArrayList<>();
        for (int off = 0; off < data.length; off += CHUNK_BYTES) {
            int len = Math.min(CHUNK_BYTES, data.length - off);
            byte[] c = new byte[len];
            System.arraycopy(data, off, c, 0, len);
            chunks.add(c);
        }
        return chunks;
    }

    // ---- query helpers ----

    private java.util.Comparator<String> nullsLast() {
        return (a, b) -> a == null ? (b == null ? 0 : 1) : (b == null ? -1 : a.compareTo(b));
    }

    private java.util.Comparator<Double> nullDoubles() {
        return (a, b) -> a == null ? (b == null ? 0 : 1) : (b == null ? -1 : Double.compare(a, b));
    }

    private boolean doubleEq(Double a, Double b) {
        if (a == null || b == null)
            return a == null && b == null;
        return a.doubleValue() == b.doubleValue() || (Double.isNaN(a) && Double.isNaN(b));
    }

    private List<String> queryStrings(String sql) throws SQLException {
        List<String> out = new ArrayList<>();
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             ResultSet rs = c.createStatement().executeQuery(sql)) {
            while (rs.next())
                out.add(rs.getString(1));
        }
        return out;
    }

    private List<Double> queryDoubles(String sql) throws SQLException {
        List<Double> out = new ArrayList<>();
        try (Connection c = DriverManager.getConnection(container.getJdbcUrl(), container.getUsername(), container.getPassword());
             ResultSet rs = c.createStatement().executeQuery(sql)) {
            while (rs.next()) {
                double d = rs.getDouble(1);
                out.add(rs.wasNull() ? null : d);
            }
        }
        return out;
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
