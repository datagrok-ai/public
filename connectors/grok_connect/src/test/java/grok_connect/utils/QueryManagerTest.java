package grok_connect.utils;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.FuncCall;
import grok_connect.resultset.ResultSetManager;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import serialization.BigIntColumn;
import serialization.DataFrame;
import serialization.Types;
import java.lang.reflect.Field;
import java.lang.reflect.Proxy;
import java.sql.ResultSetMetaData;

// int8AsInt32 option parsing and the sticky ALLOW_COL_TYPE_CHANGE tag, without a live database.
public class QueryManagerTest {
    private static final String CALL = "{\"id\":\"q\",\"func\":{\"#type\":\"DataQuery\",\"query\":\"select 1\","
            + "\"connection\":{\"dataSource\":\"Postgres\"}},\"options\":%s}";

    @BeforeAll
    public static void initProviders() {
        if (GrokConnect.providerManager == null)
            GrokConnect.providerManager = new ProviderManager();
    }

    private static FuncCall parse(String options) {
        return GrokConnect.gson.fromJson(String.format(CALL, options), FuncCall.class);
    }

    private static QueryManager manager(String options) {
        return new QueryManager(String.format(CALL, options));
    }

    private static ResultSetMetaData int8Meta() {
        return (ResultSetMetaData) Proxy.newProxyInstance(QueryManagerTest.class.getClassLoader(),
                new Class<?>[] {ResultSetMetaData.class}, (proxy, method, args) -> {
                    switch (method.getName()) {
                        case "getColumnCount": return 1;
                        case "getColumnType": return java.sql.Types.BIGINT;
                        case "getColumnTypeName":
                        case "getColumnLabel": return "int8";
                        case "getPrecision":
                        case "getScale": return 0;
                        default: throw new UnsupportedOperationException(method.getName());
                    }
                });
    }

    private static String int8TypeFor(String options) throws Exception {
        Field field = QueryManager.class.getDeclaredField("resultSetManager");
        field.setAccessible(true);
        ResultSetManager resultSetManager = (ResultSetManager) field.get(manager(options));
        resultSetManager.init(int8Meta(), 2);
        resultSetManager.processValue(10L, 1);
        return resultSetManager.getProcessedColumns()[0].getType();
    }

    private static DataFrame chunk(boolean downcast) {
        BigIntColumn col = new BigIntColumn("int8", 2);
        col.setDowncastAllowed(downcast);
        col.add(1L);
        DataFrame df = new DataFrame();
        df.addColumn(col);
        return df;
    }

    @Test
    public void int8AsInt32ReadsGsonBooleanAndString() {
        Assertions.assertTrue(QueryManager.int8AsInt32(parse("{\"int8AsInt32\":true}").options));
        Assertions.assertTrue(QueryManager.int8AsInt32(parse("{\"int8AsInt32\":\"true\"}").options));
        Assertions.assertFalse(QueryManager.int8AsInt32(parse("{}").options));
        Assertions.assertFalse(QueryManager.int8AsInt32(parse("{\"int8AsInt32\":false}").options));
        Assertions.assertFalse(QueryManager.int8AsInt32(parse("{\"int8AsInt32\":null}").options));
    }

    @Test
    public void initParamsWiresInt8AsInt32IntoTheResultSetManager() throws Exception {
        Assertions.assertEquals(Types.INT, int8TypeFor("{\"int8AsInt32\":true}"));
        Assertions.assertEquals(Types.INT, int8TypeFor("{\"int8AsInt32\":\"true\"}"));
        Assertions.assertEquals(Types.BIG_INT, int8TypeFor("{}"));
    }

    @Test
    public void fetchSizeClampsAt500000Rows() {
        QueryManager m = manager("{}");
        m.reportSerialized(100, 100);
        Assertions.assertEquals(500000, m.getFetchSize(chunk(false), m.getWireBytesPerRow()), "1 B/row -> 10 MB target = 10M rows, clamped");
        m.reportSerialized(100, 1);
        Assertions.assertEquals(100000, m.getFetchSize(chunk(false), m.getWireBytesPerRow()), "100 B/row -> 100K rows, not clamped");
        Assertions.assertEquals(100000, m.getFetchSize(chunk(false), 100f), "the caller's snapshot wins over the field");
    }

    @Test
    public void allowColTypeChangeIsDecidedOnChunkOneAndSticky() {
        QueryManager downcasting = manager("{\"int8AsInt32\":true}");
        DataFrame first = chunk(true);
        downcasting.tagChunk(first, 1);
        Assertions.assertEquals("1", first.getTags().get(QueryManager.CHUNK_NUMBER_TAG));
        Assertions.assertEquals("true", first.getTags().get(QueryManager.ALLOW_COL_TYPE_CHANGE_TAG));
        for (int n = 2; n <= 3; n++) {
            DataFrame next = chunk(false);
            downcasting.tagChunk(next, n);
            Assertions.assertEquals(String.valueOf(n), next.getTags().get(QueryManager.CHUNK_NUMBER_TAG));
            Assertions.assertEquals("true", next.getTags().get(QueryManager.ALLOW_COL_TYPE_CHANGE_TAG));
        }

        QueryManager plain = manager("{}");
        DataFrame noDowncast = chunk(false);
        plain.tagChunk(noDowncast, 1);
        Assertions.assertNull(noDowncast.getTags().get(QueryManager.ALLOW_COL_TYPE_CHANGE_TAG));
        DataFrame late = chunk(true);
        plain.tagChunk(late, 2);
        Assertions.assertNull(late.getTags().get(QueryManager.ALLOW_COL_TYPE_CHANGE_TAG));
    }
}
