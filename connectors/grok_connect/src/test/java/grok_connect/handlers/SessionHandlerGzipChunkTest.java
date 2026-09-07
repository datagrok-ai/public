package grok_connect.handlers;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.FuncCall;
import grok_connect.managers.ColumnManager;
import grok_connect.managers.datetime_column.SQLiteDateTimeColumnManager;
import grok_connect.managers.string_column.SQLiteStringColumnManager;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.providers.SQLiteDataProvider;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryCancelledByUser;
import grok_connect.utils.QueryManager;
import grok_connect.utils.QueryMonitor;
import grok_connect.utils.SettingsManager;
import org.eclipse.jetty.websocket.api.RemoteEndpoint;
import org.eclipse.jetty.websocket.api.Session;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import serialization.DataFrame;
import serialization.Types;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CompletionException;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

// Drives the real SessionHandler WS state machine against a file-backed SQLite database:
// every chunk's announce token must match what is actually sent, gzipped chunks must
// restore to the same d42 rows the plain run produces, and the two-stage fetch/serialize
// pipeline must produce exactly the bytes of the single-task sequencing.
public class SessionHandlerGzipChunkTest {
    private static final int ROWS = 250;
    private static final int BIG_ROWS = 5000;
    private static final int WIDE_ROWS = 2000;
    private static final String SQL = "SELECT id, name FROM t ORDER BY id";
    private static final String BIG_SQL = "SELECT id, name FROM big ORDER BY id";
    private static final String WIDE_SQL = "SELECT id, payload FROM wide ORDER BY id";
    private static final String NOT_OK = "NOT OK RESPONSE";
    private static final String GATED = "SQLiteGated";
    private static Path db;
    private static ProviderManager previousProviderManager;
    // While set, every column read of the gated provider parks on it (and trips `blocked` first).
    private static volatile CountDownLatch gate;
    private static volatile CountDownLatch blocked;

    private static final class GatedProvider extends SQLiteDataProvider {
        @Override
        public ResultSetManager getResultSetManager() {
            Map<String, ColumnManager<?>> managers = DefaultResultSetManager.getDefaultManagersMap();
            managers.put(Types.DATE_TIME, new SQLiteDateTimeColumnManager());
            managers.put(Types.STRING, new SQLiteStringColumnManager());
            return new DefaultResultSetManager(managers.values()) {
                @Override
                public boolean readFast(ResultSet resultSet, int index) throws SQLException {
                    park();
                    return super.readFast(resultSet, index);
                }

                @Override
                public void processValue(Object o, int index) {
                    park();
                    super.processValue(o, index);
                }
            };
        }

        private static void park() {
            CountDownLatch g = gate;
            if (g == null)
                return;
            blocked.countDown();
            try {
                g.await();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }

    private static final class Sent {
        final String announce;
        final byte[] bytes;

        Sent(String announce, byte[] bytes) {
            this.announce = announce;
            this.bytes = bytes;
        }
    }

    private static final class Driver {
        final List<String> strings = new ArrayList<>();
        final List<byte[]> binaries = new ArrayList<>();
        final FuncCall call;
        final SessionHandler handler;

        Driver(String sql, String fetchSize, Object initFetchSize, boolean compress) {
            Session session = Mockito.mock(Session.class);
            RemoteEndpoint remote = Mockito.mock(RemoteEndpoint.class);
            Mockito.when(session.getRemote()).thenReturn(remote);
            Mockito.when(remote.sendStringByFuture(Mockito.anyString())).thenAnswer(invocation -> {
                strings.add(invocation.getArgument(0));
                return null;
            });
            Mockito.when(remote.sendBytesByFuture(Mockito.any(ByteBuffer.class))).thenAnswer(invocation -> {
                ByteBuffer buffer = invocation.getArgument(0);
                byte[] copy = new byte[buffer.remaining()];
                buffer.get(copy);
                binaries.add(copy);
                return null;
            });
            call = FuncCallBuilder.fromQuery(sql);
            call.func.type = "DataQuery";
            call.func.connection = new DataConnection();
            call.func.connection.dataSource = "SQLite";
            call.func.connection.parameters.put("db", db.toString());
            call.options.put("connectFetchSize", fetchSize);
            if (initFetchSize != null)
                call.options.put("initConnectFetchSize", initFetchSize);
            if (compress)
                call.options.put("compressChunks", true);
            handler = new SessionHandler(session, true);
        }

        void query() throws Throwable {
            handler.onMessage("QUERY " + GrokConnect.gson.toJson(call));
        }

        String lastString() {
            return strings.get(strings.size() - 1);
        }

        byte[] lastBinary() {
            return binaries.get(binaries.size() - 1);
        }

        @SuppressWarnings("unchecked")
        <T> T field(String name) throws Exception {
            Field field = SessionHandler.class.getDeclaredField(name);
            field.setAccessible(true);
            return (T) field.get(handler);
        }
    }

    @BeforeAll
    @SuppressWarnings("unchecked")
    public static void init() throws Exception {
        SettingsManager.getInstance();
        previousProviderManager = GrokConnect.providerManager;
        GrokConnect.providerManager = new ProviderManager();
        // SQLite is not in ProviderManager.PROVIDER_CLASSES (file-based, never advertised), but its
        // driver ships with grok_connect — the only database this test can reach without Docker.
        Field providersField = ProviderManager.class.getDeclaredField("providersMap");
        providersField.setAccessible(true);
        Map<String, JdbcDataProvider> providers = (Map<String, JdbcDataProvider>) providersField.get(GrokConnect.providerManager);
        providers.put("SQLite", new SQLiteDataProvider());
        providers.put(GATED, new GatedProvider());
        Field poolField = GrokConnect.class.getDeclaredField("threadPool");
        poolField.setAccessible(true);
        if (poolField.get(null) == null)
            poolField.set(null, Executors.newCachedThreadPool());
        db = Files.createTempFile("gzip-chunks", ".sqlite");
        Class.forName("org.sqlite.JDBC");
        try (Connection c = DriverManager.getConnection("jdbc:sqlite:" + db);
             Statement s = c.createStatement()) {
            s.execute("CREATE TABLE t (id int, name text)");
            for (int i = 1; i <= ROWS; i++)
                s.execute(String.format("INSERT INTO t VALUES (%d, 'name-%d')", i, i % 17));
            s.execute("CREATE TABLE big (id int, name text)");
            c.setAutoCommit(false);
            for (int i = 1; i <= BIG_ROWS; i++)
                s.execute(String.format("INSERT INTO big VALUES (%d, 'name-%d')", i, i % 101));
            // Distinct fixed-width ~5 KB rows: serialized bytes/row is the same for every chunk, so auto sizing
            // must land on the same row counts whichever chunk's measurement feeds it.
            s.execute("CREATE TABLE wide (id int, payload text)");
            for (int i = 1; i <= WIDE_ROWS; i++)
                s.execute(String.format("INSERT INTO wide VALUES (%d, hex(randomblob(2500)))", i));
            c.commit();
        }
    }

    @AfterAll
    public static void cleanup() throws IOException {
        GrokConnect.providerManager = previousProviderManager;
        Files.deleteIfExists(db);
    }

    private static byte[] gunzip(byte[] data) throws IOException {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        try (GZIPInputStream in = new GZIPInputStream(new ByteArrayInputStream(data))) {
            byte[] buf = new byte[8192];
            int n;
            while ((n = in.read(buf)) > 0)
                out.write(buf, 0, n);
        }
        return out.toByteArray();
    }

    private static List<Sent> drive(Object initFetchSize, boolean compress) throws Throwable {
        return drive(initFetchSize, compress, -1);
    }

    /** Runs the scenario with the fetch pipeline off and on; both must announce and send identical bytes. */
    private static List<Sent> drive(Object initFetchSize, boolean compress, int nackChunk) throws Throwable {
        return bothModes(SQL, "100", initFetchSize, compress, nackChunk);
    }

    private static List<Sent> bothModes(String sql, String fetchSize, Object initFetchSize, boolean compress, int nackChunk) throws Throwable {
        List<Sent> sequential = drive(sql, fetchSize, initFetchSize, compress, nackChunk, false);
        List<Sent> pipelined = drive(sql, fetchSize, initFetchSize, compress, nackChunk, true);
        Assertions.assertEquals(sequential.size(), pipelined.size());
        for (int i = 0; i < sequential.size(); i++) {
            Assertions.assertEquals(sequential.get(i).announce, pipelined.get(i).announce, "announce of chunk " + (i + 1));
            Assertions.assertArrayEquals(sequential.get(i).bytes, pipelined.get(i).bytes, "bytes of chunk " + (i + 1));
        }
        return pipelined;
    }

    /** Plays the Datlas half of the protocol: QUERY, then SIZE RECEIVED / PART OK per chunk until COMPLETED;
     *  chunk {@code nackChunk} (0-based) is first answered with NOT OK so that GC re-announces and resends it. */
    private static List<Sent> drive(String sql, String fetchSize, Object initFetchSize, boolean compress, int nackChunk,
                                    boolean pipeline) throws Throwable {
        boolean previous = SessionHandler.PIPELINE;
        SessionHandler.PIPELINE = pipeline;
        Driver d = new Driver(sql, fetchSize, initFetchSize, compress);
        try {
            d.query();
            List<Sent> sent = new ArrayList<>();
            while (!d.lastString().startsWith("COMPLETED ")) {
                String announce = d.lastString();
                d.handler.onMessage("DATAFRAME PART SIZE RECEIVED");
                byte[] bytes = d.lastBinary();
                if (sent.size() == nackChunk) {
                    d.handler.onMessage(NOT_OK);
                    Assertions.assertEquals(announce, d.lastString());
                    d.handler.onMessage("DATAFRAME PART SIZE RECEIVED");
                    Assertions.assertArrayEquals(bytes, d.lastBinary());
                }
                sent.add(new Sent(announce, bytes));
                d.handler.onMessage("PART OK");
            }
            Assertions.assertEquals("COMPLETED " + sent.size(), d.lastString());
            Assertions.assertEquals(sent.size() + (nackChunk >= 0 ? 1 : 0), d.binaries.size());
            d.handler.onMessage("COMPLETED_OK");
            return sent;
        } finally {
            d.handler.onClose();
            SessionHandler.PIPELINE = previous;
        }
    }

    private static int announcedSize(String announce) {
        return Integer.parseInt(announce.substring("DATAFRAME PART SIZE: ".length()).split(" ")[0]);
    }

    private static DataFrame frameOf(Sent sent) throws IOException {
        boolean gzipped = sent.announce.endsWith(" gzip=true");
        return DataFrame.fromByteArray(gzipped ? gunzip(sent.bytes) : sent.bytes);
    }

    private static int rowsOf(Sent sent) throws IOException {
        return frameOf(sent).rowCount;
    }

    @Test
    public void withoutTheOptionChunksArePlainAndUntagged() throws Throwable {
        List<Sent> sent = drive("100000", false);
        Assertions.assertEquals(1, sent.size());
        Assertions.assertEquals("DATAFRAME PART SIZE: " + sent.get(0).bytes.length, sent.get(0).announce);
        Assertions.assertEquals(ROWS, DataFrame.fromByteArray(sent.get(0).bytes).rowCount);
    }

    @Test
    public void everyChunkGzippedWhenInitFetchSizeIsAString() throws Throwable {
        List<Sent> sent = drive("100", true);
        Assertions.assertEquals(3, sent.size());
        int rows = 0;
        for (Sent chunk : sent) {
            Assertions.assertEquals("DATAFRAME PART SIZE: " + chunk.bytes.length + " gzip=true", chunk.announce);
            Assertions.assertEquals(chunk.bytes.length, announcedSize(chunk.announce));
            rows += rowsOf(chunk);
        }
        Assertions.assertEquals(ROWS, rows);
    }

    @Test
    public void chunkOneStaysPlainWhenInitFetchSizeIs100() throws Throwable {
        List<Sent> sent = drive(100, true);
        Assertions.assertEquals(3, sent.size());
        Assertions.assertEquals("DATAFRAME PART SIZE: " + sent.get(0).bytes.length, sent.get(0).announce);
        Assertions.assertEquals(100, DataFrame.fromByteArray(sent.get(0).bytes).rowCount);
        Assertions.assertTrue(sent.get(1).announce.endsWith(" gzip=true"));
        Assertions.assertTrue(sent.get(2).announce.endsWith(" gzip=true"));
        Assertions.assertEquals(ROWS, rowsOf(sent.get(0)) + rowsOf(sent.get(1)) + rowsOf(sent.get(2)));
    }

    @Test
    public void nackedChunkIsReannouncedAndResentUnchanged() throws Throwable {
        List<Sent> sent = drive("100", true, 1);
        Assertions.assertEquals(3, sent.size());
        Assertions.assertTrue(sent.get(1).announce.endsWith(" gzip=true"));
        Assertions.assertEquals(sent.get(1).bytes.length, announcedSize(sent.get(1).announce));
        Assertions.assertEquals(ROWS, rowsOf(sent.get(0)) + rowsOf(sent.get(1)) + rowsOf(sent.get(2)));
    }

    @Test
    public void chunkOneStaysPlainWhenInitFetchSizeIsAbsent() throws Throwable {
        List<Sent> sent = drive(null, true);
        Assertions.assertEquals(3, sent.size());
        Assertions.assertFalse(sent.get(0).announce.contains("gzip="));
        Assertions.assertTrue(sent.get(1).announce.endsWith(" gzip=true"));
        Assertions.assertEquals(ROWS, rowsOf(sent.get(0)) + rowsOf(sent.get(1)) + rowsOf(sent.get(2)));
    }

    @Test
    public void tenChunksAreByteIdenticalAndOrderedInBothModes() throws Throwable {
        List<Sent> sent = bothModes(BIG_SQL, "500", "500", true, -1);
        Assertions.assertEquals(10, sent.size());
        int nextId = 1;
        for (Sent chunk : sent) {
            DataFrame df = frameOf(chunk);
            Assertions.assertEquals(500, df.rowCount);
            Assertions.assertEquals(String.valueOf(sent.indexOf(chunk) + 1), df.getTag(QueryManager.CHUNK_NUMBER_TAG));
            for (int i = 0; i < df.rowCount; i++)
                Assertions.assertEquals(nextId++, ((Number) df.getColumn(0).get(i)).intValue());
        }
        Assertions.assertEquals(BIG_ROWS + 1, nextId);
    }

    @Test
    public void cancelledFetchPropagatesThroughThePipeline() throws Throwable {
        boolean previous = SessionHandler.PIPELINE;
        SessionHandler.PIPELINE = true;
        Driver d = new Driver(SQL, "100", "100", true);
        try {
            d.query();
            d.<CompletableFuture<DataFrame>>field("fetchFuture").join();
            QueryMonitor.getInstance().addCancelledResultSet(d.call.id);
            d.handler.onMessage("DATAFRAME PART SIZE RECEIVED");
            d.handler.onMessage("PART OK");
            Assertions.assertTrue(d.lastString().startsWith("DATAFRAME PART SIZE: "));
            d.handler.onMessage("DATAFRAME PART SIZE RECEIVED");
            Assertions.assertThrows(QueryCancelledByUser.class, () -> d.handler.onMessage("PART OK"));
            Assertions.assertEquals(2, d.binaries.size());
        } finally {
            QueryMonitor.getInstance().removeResultSet(d.call.id);
            d.handler.onClose();
            SessionHandler.PIPELINE = previous;
        }
    }

    @Test
    public void autoSizedChunksHaveIdenticalRowCountsInBothModes() throws Throwable {
        List<Sent> sequential = drive(WIDE_SQL, "1 MB", "150", false, -1, false);
        List<Sent> pipelined = drive(WIDE_SQL, "1 MB", "150", false, -1, true);
        Assertions.assertEquals(sequential.size(), pipelined.size());
        Assertions.assertTrue(sequential.size() >= 4, "auto sizing must split the table, got " + sequential.size() + " chunks");
        int rows = 0;
        for (int i = 0; i < sequential.size(); i++) {
            Assertions.assertEquals(rowsOf(sequential.get(i)), rowsOf(pipelined.get(i)), "rows of chunk " + (i + 1));
            rows += rowsOf(pipelined.get(i));
        }
        Assertions.assertEquals(WIDE_ROWS, rows);
        Assertions.assertEquals(150, rowsOf(pipelined.get(0)));
        Assertions.assertNotEquals(150, rowsOf(pipelined.get(2)), "chunk 3 must be sized from serialized bytes, not the init size");
    }

    @Test
    public void closeMidStreamCancelsBothStagesAndReleasesTheConnection() throws Throwable {
        boolean previous = SessionHandler.PIPELINE;
        SessionHandler.PIPELINE = true;
        Driver d = new Driver(BIG_SQL, "500", "500", true);
        d.call.func.connection.dataSource = GATED;
        CountDownLatch release = new CountDownLatch(1);
        try {
            d.query();
            d.<CompletableFuture<DataFrame>>field("fetchFuture").join();
            blocked = new CountDownLatch(1);
            gate = release;
            d.handler.onMessage("DATAFRAME PART SIZE RECEIVED");
            Assertions.assertTrue(blocked.await(10, TimeUnit.SECONDS), "fetch(3) must be parked mid-row");
            CompletableFuture<?> serializing = d.field("completableFuture");
            CompletableFuture<?> fetching = d.field("fetchFuture");
            Assertions.assertFalse(fetching.isDone());
            long releaseDelayMs = 300;
            new Thread(() -> {
                try {
                    Thread.sleep(releaseDelayMs);
                } catch (InterruptedException ignore) {}
                release.countDown();
            }).start();
            long start = System.nanoTime();
            d.handler.onClose();
            Assertions.assertTrue((System.nanoTime() - start) / 1_000_000 >= releaseDelayMs - 50, "onClose must wait for the parked row");
            Assertions.assertTrue(serializing.isDone());
            Assertions.assertTrue(fetching.isCompletedExceptionally());
            Throwable cause = Assertions.assertThrows(CompletionException.class, fetching::join);
            while (cause.getCause() != null)
                cause = cause.getCause();
            Assertions.assertTrue(cause instanceof QueryCancelledByUser, "fetch must end in QueryCancelledByUser, got " + cause);
            QueryManager queryManager = d.field("queryManager");
            Field connection = QueryManager.class.getDeclaredField("connection");
            connection.setAccessible(true);
            Assertions.assertTrue(((Connection) connection.get(queryManager)).isClosed());
            Assertions.assertFalse(QueryMonitor.getInstance().checkCancelledIdResultSet(d.call.id));
        } finally {
            gate = null;
            release.countDown();
            QueryMonitor.getInstance().removeResultSet(d.call.id);
            SessionHandler.PIPELINE = previous;
        }
    }
}
