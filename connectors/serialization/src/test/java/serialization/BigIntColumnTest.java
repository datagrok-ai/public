package serialization;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

// BigIntColumn dual storage: int-mode downcast (opt-in, sticky, lossless) vs today's
// string-dictionary bigint payload.
public class BigIntColumnTest {
    private static BigIntColumn column(boolean downcast, Object... values) {
        BigIntColumn col = new BigIntColumn("c", 2);
        col.setDowncastAllowed(downcast);
        for (Object value : values)
            col.add(value);
        return col;
    }

    private static byte[] encode(Column<?> col) {
        BufferAccessor buf = new BufferAccessor();
        col.encode(buf);
        return buf.toUint8List();
    }

    @Test
    public void intModeWritesExactlyIntColumnPayload() {
        BigIntColumn col = column(true, 1L, null, -5L, (long) Integer.MAX_VALUE, (long) Integer.MIN_VALUE + 1, 7);
        Assertions.assertEquals(Types.INT, col.getType());
        Assertions.assertEquals(6, col.getLength());
        Assertions.assertEquals(1L, col.get(0));
        Assertions.assertNull(col.get(1));
        Assertions.assertTrue(col.isNone(1));
        Assertions.assertEquals(7L, col.get(5));
        IntColumn ints = new IntColumn("c", new Integer[] {1, null, -5, Integer.MAX_VALUE, Integer.MIN_VALUE + 1, 7});
        Assertions.assertArrayEquals(encode(ints), encode(col));
    }

    @Test
    public void bigintModeWritesTodaysStringPayload() {
        BigIntColumn col = column(false, 1L, null, -5L, 1L << 40);
        Assertions.assertEquals(Types.BIG_INT, col.getType());
        BigIntColumn strings = new BigIntColumn("c", new String[] {"1", "", "-5", String.valueOf(1L << 40)});
        Assertions.assertEquals(Types.BIG_INT, strings.getType());
        Assertions.assertArrayEquals(encode(strings), encode(col));
    }

    @Test
    public void stickySetExplicitlyReportsAndEncodesBigint() {
        BigIntColumn col = new BigIntColumn("c", 2);
        col.setDowncastAllowed(true);
        col.setSticky(true);
        col.add(1L);
        Assertions.assertEquals(Types.BIG_INT, col.getType());
        Assertions.assertArrayEquals(encode(new BigIntColumn("c", new String[] {"1"})), encode(col));
    }

    @Test
    public void downcastNotAllowedStaysBigint() {
        Assertions.assertEquals(Types.BIG_INT, column(false, 1L, 2L).getType());
    }

    @Test
    public void overflowIsStickyAcrossEmpty() {
        BigIntColumn col = column(true, 1L, 2L);
        Assertions.assertEquals(Types.INT, col.getType());
        col.add(1L << 40);
        Assertions.assertEquals(Types.BIG_INT, col.getType());
        Assertions.assertTrue(col.isSticky());
        col.empty();
        Assertions.assertEquals(0, col.getLength());
        col.add(3L);
        Assertions.assertEquals(Types.BIG_INT, col.getType());
        Assertions.assertEquals(3L, col.get(0));
    }

    @Test
    public void intMinValueIsNotDowncast() {
        Assertions.assertEquals(Types.BIG_INT, column(true, (long) Integer.MIN_VALUE).getType());
        Assertions.assertEquals(Types.BIG_INT, column(true, (long) Integer.MAX_VALUE + 1).getType());
    }

    @Test
    public void emptyAndAllNoneColumns() {
        BigIntColumn empty = column(true);
        Assertions.assertEquals(Types.BIG_INT, empty.getType());
        BigIntColumn nones = column(true, null, null);
        Assertions.assertEquals(Types.BIG_INT, nones.getType());
        Assertions.assertTrue(nones.isNone(0) && nones.isNone(1));
    }

    @Test
    public void nonIntegralValueSwitchesToStrings() {
        BigIntColumn col = column(true, 5L, null);
        col.add("123456789012345678901234567890");
        col.add(null);
        Assertions.assertEquals(Types.BIG_INT, col.getType());
        Assertions.assertEquals(4, col.getLength());
        Assertions.assertEquals("5", col.get(0));
        Assertions.assertTrue(col.isNone(1));
        Assertions.assertEquals("123456789012345678901234567890", col.get(2));
        Assertions.assertTrue(col.isNone(3));
        col.empty();
        col.add(1L);
        Assertions.assertEquals(Types.BIG_INT, col.getType());
    }

    @Test
    public void longMinValueIsAValueNotNone() {
        BigIntColumn col = column(true, Long.MIN_VALUE);
        Assertions.assertFalse(col.isNone(0));
        Assertions.assertEquals("-9223372036854775808", col.get(0));
        Assertions.assertEquals(Types.BIG_INT, col.getType());
    }

    @Test
    public void intModeRoundTripsAsIntColumn() {
        DataFrame df = new DataFrame();
        df.addColumn(column(true, 1L, null, -5L));
        DataFrame decoded = DataFrame.fromByteArray(df.toByteArray());
        Column<?> col = decoded.getColumn("c");
        Assertions.assertTrue(col instanceof IntColumn);
        Assertions.assertEquals(1, col.get(0));
        Assertions.assertTrue(col.isNone(1));
        Assertions.assertEquals(-5, col.get(2));
    }

    @Test
    public void bigintModeRoundTripsAsStrings() {
        DataFrame df = new DataFrame();
        df.addColumn(column(true, 1L, null, 1L << 40));
        DataFrame decoded = DataFrame.fromByteArray(df.toByteArray());
        BigIntColumn col = (BigIntColumn) decoded.getColumn("c");
        Assertions.assertEquals(Types.BIG_INT, col.getType());
        Assertions.assertEquals("1", col.get(0));
        Assertions.assertTrue(col.isNone(1));
        Assertions.assertEquals(String.valueOf(1L << 40), col.get(2));
    }
}
