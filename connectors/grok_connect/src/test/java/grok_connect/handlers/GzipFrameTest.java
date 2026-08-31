package grok_connect.handlers;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryManager;
import org.eclipse.jetty.websocket.api.Session;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.Map;
import java.util.Random;
import java.util.zip.GZIPInputStream;

// Gzip-at-GC building blocks without a database: the gzip helper, the option parsing,
// the chunk-1 rule, and serialize() producing a Frame the standard GZIPInputStream restores.
public class GzipFrameTest {
    private static final String CALL = "{\"id\":\"q\",\"func\":{\"#type\":\"DataQuery\",\"query\":\"select 1\","
            + "\"connection\":{\"dataSource\":\"Postgres\"}},\"options\":%s}";

    @BeforeAll
    public static void initProviders() {
        if (GrokConnect.providerManager == null)
            GrokConnect.providerManager = new ProviderManager();
    }

    private static Map<String, Object> options(String json) {
        return GrokConnect.gson.fromJson(String.format(CALL, json), FuncCall.class).options;
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

    private static DataFrame frame(int rows) {
        IntColumn ints = new IntColumn("i", rows);
        StringColumn strings = new StringColumn("s", rows);
        for (int i = 0; i < rows; i++) {
            ints.add(i * 7);
            strings.add("value-" + (i % 13));
        }
        DataFrame df = new DataFrame();
        df.addColumn(ints);
        df.addColumn(strings);
        return df;
    }

    private static SessionHandler handler(String options) throws Exception {
        SessionHandler handler = new SessionHandler(Mockito.mock(Session.class), true);
        Field field = SessionHandler.class.getDeclaredField("queryManager");
        field.setAccessible(true);
        field.set(handler, new QueryManager(String.format(CALL, options)));
        return handler;
    }

    @Test
    public void gzipRoundTripsAtEveryLevel() throws IOException {
        byte[] data = new byte[100_000];
        new Random(42).nextBytes(data);
        for (int i = 0; i < data.length; i += 2)
            data[i] = 0;
        for (int level = 1; level <= 9; level++)
            Assertions.assertArrayEquals(data, gunzip(SessionHandler.gzip(data, level)), "level " + level);
        Assertions.assertTrue(SessionHandler.gzip(data, 9).length <= SessionHandler.gzip(data, 1).length);
    }

    @Test
    public void gzipLevelDefaultsToOneAndClampsToOneToNine() {
        Assertions.assertEquals(1, QueryManager.gzipLevel(options("{}")));
        Assertions.assertEquals(6, QueryManager.gzipLevel(options("{\"gzipLevel\":6}")));
        Assertions.assertEquals(3, QueryManager.gzipLevel(options("{\"gzipLevel\":\"3\"}")));
        Assertions.assertEquals(1, QueryManager.gzipLevel(options("{\"gzipLevel\":0}")));
        Assertions.assertEquals(1, QueryManager.gzipLevel(options("{\"gzipLevel\":-5}")));
        Assertions.assertEquals(9, QueryManager.gzipLevel(options("{\"gzipLevel\":42}")));
    }

    @Test
    public void compressChunksReadsGsonBooleanAndString() {
        Assertions.assertTrue(QueryManager.compressChunks(options("{\"compressChunks\":true}")));
        Assertions.assertTrue(QueryManager.compressChunks(options("{\"compressChunks\":\"true\"}")));
        Assertions.assertFalse(QueryManager.compressChunks(options("{}")));
        Assertions.assertFalse(QueryManager.compressChunks(options("{\"compressChunks\":false}")));
        Assertions.assertFalse(QueryManager.compressChunks(options("{\"compressChunks\":null}")));
    }

    @Test
    public void chunkOneIsGzippedOnlyWhenInitFetchSizeGivenAndNot100() {
        Assertions.assertFalse(QueryManager.gzipChunk(options("{\"compressChunks\":true}"), 1));
        Assertions.assertFalse(QueryManager.gzipChunk(options("{\"compressChunks\":true,\"initConnectFetchSize\":100}"), 1));
        Assertions.assertTrue(QueryManager.gzipChunk(options("{\"compressChunks\":true,\"initConnectFetchSize\":\"100000\"}"), 1));
        Assertions.assertTrue(QueryManager.gzipChunk(options("{\"compressChunks\":true,\"initConnectFetchSize\":\"100\"}"), 1));
        Assertions.assertTrue(QueryManager.gzipChunk(options("{\"compressChunks\":true,\"initConnectFetchSize\":200}"), 1));
        for (int n = 2; n <= 3; n++) {
            Assertions.assertTrue(QueryManager.gzipChunk(options("{\"compressChunks\":true}"), n));
            Assertions.assertTrue(QueryManager.gzipChunk(options("{\"compressChunks\":true,\"initConnectFetchSize\":100}"), n));
        }
        Assertions.assertFalse(QueryManager.gzipChunk(options("{}"), 2));
        Assertions.assertFalse(QueryManager.gzipChunk(options("{\"compressChunks\":false,\"initConnectFetchSize\":\"100000\"}"), 1));
    }

    @Test
    public void serializeProducesPlainFrameWithoutTheOption() throws Exception {
        DataFrame df = frame(500);
        SessionHandler.Frame frame = handler("{}").serialize(df, 2);
        Assertions.assertFalse(frame.gzipped);
        Assertions.assertEquals(500, frame.rowCount);
        Assertions.assertArrayEquals(df.toByteArray(), frame.bytes);
    }

    @Test
    public void serializeGzipsTheD42BytesWhenOptedIn() throws Exception {
        DataFrame df = frame(500);
        byte[] plain = df.toByteArray();
        SessionHandler.Frame frame = handler("{\"compressChunks\":true,\"gzipLevel\":1}").serialize(df, 2);
        Assertions.assertTrue(frame.gzipped);
        Assertions.assertEquals(500, frame.rowCount);
        Assertions.assertTrue(frame.bytes.length < plain.length);
        Assertions.assertArrayEquals(plain, gunzip(frame.bytes));
        Assertions.assertEquals(500, DataFrame.fromByteArray(gunzip(frame.bytes)).rowCount);
    }

    @Test
    public void serializeAppliesTheChunkOneRule() throws Exception {
        SessionHandler handler = handler("{\"compressChunks\":true,\"initConnectFetchSize\":100}");
        Assertions.assertFalse(handler.serialize(frame(10), 1).gzipped);
        Assertions.assertTrue(handler.serialize(frame(10), 2).gzipped);
        Assertions.assertTrue(handler("{\"compressChunks\":true,\"initConnectFetchSize\":\"100000\"}").serialize(frame(10), 1).gzipped);
    }
}
