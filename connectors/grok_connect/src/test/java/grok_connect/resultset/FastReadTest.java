package grok_connect.resultset;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import serialization.BigIntColumn;
import serialization.BoolColumn;
import serialization.Column;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import serialization.Types;
import java.lang.reflect.Proxy;
import java.math.BigDecimal;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.Timestamp;

// Typed-accessor row loop against an in-memory JDBC ResultSet: the fast path must keep
// getObject's null and narrowing semantics exactly.
public class FastReadTest {
    private static final String[] NAMES = {"int4", "float8", "float4", "numeric", "bool", "text", "int8", "ts"};
    private static final int[] TYPES = {java.sql.Types.INTEGER, java.sql.Types.DOUBLE, java.sql.Types.REAL,
            java.sql.Types.NUMERIC, java.sql.Types.BOOLEAN, java.sql.Types.VARCHAR, java.sql.Types.BIGINT, java.sql.Types.TIMESTAMP};
    private static final Object[][] ROWS = {
            {1, 1.5, 2.5f, new BigDecimal("3.250"), true, "a", 10L, new Timestamp(0)},
            {null, null, null, null, null, null, null, null},
            {-7, 1e300, -0.5f, new BigDecimal("-1"), false, "", 20L, new Timestamp(1000)},
    };

    private static ResultSetMetaData meta() {
        return (ResultSetMetaData) Proxy.newProxyInstance(FastReadTest.class.getClassLoader(),
                new Class<?>[] {ResultSetMetaData.class}, (proxy, method, args) -> {
                    switch (method.getName()) {
                        case "getColumnCount": return NAMES.length;
                        case "getColumnType": return TYPES[(Integer) args[0] - 1];
                        case "getColumnTypeName":
                        case "getColumnLabel": return NAMES[(Integer) args[0] - 1];
                        case "getPrecision":
                        case "getScale": return 0;
                        default: throw new UnsupportedOperationException(method.getName());
                    }
                });
    }

    private static ResultSet resultSet(Object[][] rows) {
        int[] row = {-1};
        boolean[] wasNull = {false};
        return (ResultSet) Proxy.newProxyInstance(FastReadTest.class.getClassLoader(),
                new Class<?>[] {ResultSet.class}, (proxy, method, args) -> {
                    if (method.getName().equals("next"))
                        return ++row[0] < rows.length;
                    if (method.getName().equals("wasNull"))
                        return wasNull[0];
                    if (method.getName().equals("getMetaData"))
                        return meta();
                    Object value = rows[row[0]][(Integer) args[0] - 1];
                    wasNull[0] = value == null;
                    switch (method.getName()) {
                        case "getObject": return value;
                        case "getInt": return value == null ? 0 : ((Number) value).intValue();
                        case "getDouble": return value == null ? 0.0 : ((Number) value).doubleValue();
                        case "getFloat": return value == null ? 0f : ((Number) value).floatValue();
                        case "getBoolean": return value != null && (Boolean) value;
                        case "getString": return value == null ? null : value.toString();
                        default: throw new UnsupportedOperationException(method.getName());
                    }
                });
    }

    private static boolean[] fill(ResultSetManager manager, Object[][] rows) throws Exception {
        ResultSet rs = resultSet(rows);
        boolean[] fast = new boolean[NAMES.length];
        while (rs.next())
            for (int c = 1; c <= NAMES.length; c++) {
                fast[c - 1] = manager.readFast(rs, c);
                if (!fast[c - 1])
                    manager.processValue(rs.getObject(c), c);
            }
        return fast;
    }

    @Test
    public void fastPathKeepsGetObjectSemantics() throws Exception {
        ResultSetManager manager = DefaultResultSetManager.getDefaultManager();
        manager.setInt8AsInt32(true);
        manager.init(meta(), 2);
        boolean[] fast = fill(manager, ROWS);
        Assertions.assertArrayEquals(new boolean[] {true, true, true, false, true, true, false, false}, fast);
        Column<?>[] cols = manager.getProcessedColumns();

        IntColumn int4 = (IntColumn) cols[0];
        Assertions.assertEquals(1, int4.get(0));
        Assertions.assertTrue(int4.isNone(1));
        Assertions.assertEquals(-7, int4.get(2));

        FloatColumn float8 = (FloatColumn) cols[1];
        Assertions.assertEquals(1.5f, float8.get(0));
        Assertions.assertTrue(float8.isNone(1));
        Assertions.assertEquals(Float.POSITIVE_INFINITY, float8.get(2));

        FloatColumn float4 = (FloatColumn) cols[2];
        Assertions.assertEquals(2.5f, float4.get(0));
        Assertions.assertTrue(float4.isNone(1));
        Assertions.assertEquals(-0.5f, float4.get(2));

        FloatColumn numeric = (FloatColumn) cols[3];
        Assertions.assertEquals(3.25f, numeric.get(0));
        Assertions.assertTrue(numeric.isNone(1));
        Assertions.assertEquals(-1f, numeric.get(2));

        BoolColumn bool = (BoolColumn) cols[4];
        Assertions.assertTrue(bool.get(0));
        Assertions.assertFalse(bool.get(1));
        Assertions.assertFalse(bool.get(2));

        StringColumn text = (StringColumn) cols[5];
        Assertions.assertEquals("a", text.get(0));
        Assertions.assertEquals("", text.get(1));
        Assertions.assertTrue(text.isNone(1));
        Assertions.assertEquals("", text.get(2));

        BigIntColumn int8 = (BigIntColumn) cols[6];
        Assertions.assertEquals(Types.INT, int8.getType());
        Assertions.assertEquals(10L, int8.get(0));
        Assertions.assertTrue(int8.isNone(1));
        Assertions.assertEquals(20L, int8.get(2));

        Assertions.assertEquals(Types.DATE_TIME, cols[7].getType());
        Assertions.assertEquals(3, cols[7].getLength());
    }

    @Test
    public void int8DowncastIsOptInAndStickyAcrossChunks() throws Exception {
        ResultSetManager manager = DefaultResultSetManager.getDefaultManager();
        manager.init(meta(), 2);
        fill(manager, ROWS);
        Assertions.assertEquals(Types.BIG_INT, manager.getProcessedColumns()[6].getType());

        manager = DefaultResultSetManager.getDefaultManager();
        manager.setInt8AsInt32(true);
        manager.init(meta(), 2);
        fill(manager, ROWS);
        BigIntColumn int8 = (BigIntColumn) manager.getProcessedColumns()[6];
        Assertions.assertEquals(Types.INT, int8.getType());

        manager.detach(2);
        Object[][] overflow = {{1, 1.5, 2.5f, new BigDecimal("1"), true, "b", 1L << 40, new Timestamp(0)}};
        fill(manager, overflow);
        BigIntColumn overflowed = (BigIntColumn) manager.getProcessedColumns()[6];
        Assertions.assertNotSame(int8, overflowed);
        Assertions.assertEquals(Types.INT, int8.getType());
        Assertions.assertEquals(ROWS.length, int8.getLength());
        Assertions.assertEquals(Types.BIG_INT, overflowed.getType());
        Assertions.assertEquals(1L << 40, overflowed.get(0));

        manager.detach(2);
        fill(manager, ROWS);
        BigIntColumn sticky = (BigIntColumn) manager.getProcessedColumns()[6];
        Assertions.assertEquals(Types.BIG_INT, sticky.getType());
        Assertions.assertEquals(10L, sticky.get(0));

        for (int hop = 0; hop < 2; hop++) {
            manager.detach(2);
            fill(manager, ROWS);
            BigIntColumn later = (BigIntColumn) manager.getProcessedColumns()[6];
            Assertions.assertEquals(Types.BIG_INT, later.getType(), "detach hop " + (hop + 3));
            Assertions.assertTrue(later.isSticky());
        }
    }

    @Test
    public void detachLeavesFilledColumnsToTheirDataFrame() throws Exception {
        ResultSetManager manager = DefaultResultSetManager.getDefaultManager();
        manager.init(meta(), 2);
        fill(manager, ROWS);
        Column<?>[] filled = manager.getProcessedColumns().clone();
        manager.detach(5);
        Column<?>[] fresh = manager.getProcessedColumns();
        for (int i = 0; i < filled.length; i++) {
            Assertions.assertNotSame(filled[i], fresh[i]);
            Assertions.assertSame(filled[i].getClass(), fresh[i].getClass());
            Assertions.assertEquals(filled[i].getName(), fresh[i].getName());
            Assertions.assertEquals(ROWS.length, filled[i].getLength());
            Assertions.assertEquals(0, fresh[i].getLength());
            Assertions.assertEquals(5, fresh[i].getInitColumnSize());
        }
        Assertions.assertEquals(1, filled[0].get(0));
        Assertions.assertEquals("a", filled[5].get(0));
        fill(manager, ROWS);
        Assertions.assertEquals(ROWS.length, fresh[0].getLength());
        Assertions.assertEquals(ROWS.length, filled[0].getLength());
    }
}
