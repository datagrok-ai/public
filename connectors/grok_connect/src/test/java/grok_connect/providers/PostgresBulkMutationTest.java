package grok_connect.providers;

import java.nio.charset.StandardCharsets;
import java.sql.Connection;
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
import grok_connect.utils.ProviderManager;
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
        for (String table : new String[] {"mut_bulk_copy", "mut_bulk_batch", "mut_bulk_abort", "mut_bulk_bad"})
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
}
