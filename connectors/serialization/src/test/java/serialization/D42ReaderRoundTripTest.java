package serialization;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

// Java writer -> Java reader round-trip coverage for the d42 reader (WO-1).
// No Docker / no Dart: builds Java DataFrames, serializes with the existing
// writer, decodes with the new reader, and asserts every cell.
public class D42ReaderRoundTripTest {

    // Astral-plane emoji (U+1F600) - UTF-16 surrogate pair, exercises the
    // UTF-16-code-unit string-list slicing.
    private static final String EMOJI = "\uD83D\uDE00x";
    // CJK (Japanese) characters.
    private static final String CJK = "\u65E5\u672C\u8A9E";
    private static final String TRICKY = "with,comma\"quote\nnewline";

    private static final Integer[] INTS = {1, null, Integer.MIN_VALUE + 1, -5, 100, 0};
    private static final String[] STRS = {"hello", "", null, CJK, EMOJI, TRICKY};
    private static final Boolean[] BOOLS = {true, false, true, null, false, true};
    private static final Float[] FLOATS = {1.5F, null, Float.NaN, -3.25F, 0.0F, 1e30F};
    private static final Double[] DOUBLES = {1.5, null, 1.0000000000000002, -2.0, Double.NaN, 1e308};
    private static final Double[] DATES = {1000.0, null, -1_000_000_000.0, 1.5e15, 0.0, -1.0};
    private static final String[] BIGS = {"1152921504606846983", null, "-42", "", "0", "999"};

    private static DataFrame buildFullFrame() {
        DataFrame df = new DataFrame();
        df.name = "full";
        df.setTag("frameKey", "frameVal");

        IntColumn ints = new IntColumn("ints", INTS);
        ints.setTag("colKey", "colVal");
        df.addColumn(ints);
        df.addColumn(new StringColumn("strs", STRS));
        df.addColumn(new BoolColumn("bools", BOOLS));
        df.addColumn(new FloatColumn("floats", FLOATS));
        df.addColumn(FloatColumn.double64("doubles", DOUBLES));
        df.addColumn(new DateTimeColumn("dates", DATES));
        df.addColumn(new BigIntColumn("bigs", BIGS));
        return df;
    }

    private static void assertNoneAgrees(Column<?> col, int row, boolean expectedNone) {
        Assertions.assertEquals(expectedNone, col.isNone(row),
                "isNone mismatch at " + col.getName() + "[" + row + "]");
    }

    private static void assertFullFrame(DataFrame df) {
        Assertions.assertEquals("full", df.name);
        Assertions.assertEquals(Integer.valueOf(6), df.rowCount);
        Assertions.assertEquals("frameVal", df.getTag("frameKey"));

        IntColumn ints = (IntColumn) df.getColumn("ints");
        Assertions.assertEquals(Types.INT, ints.getType());
        Assertions.assertEquals("colVal", ints.getTag("colKey"));
        for (int r = 0; r < INTS.length; r++) {
            boolean none = INTS[r] == null;
            assertNoneAgrees(ints, r, none);
            if (!none)
                Assertions.assertEquals(INTS[r], ints.get(r));
        }

        StringColumn strs = (StringColumn) df.getColumn("strs");
        Assertions.assertEquals(Types.STRING, strs.getType());
        for (int r = 0; r < STRS.length; r++) {
            boolean none = STRS[r] == null || STRS[r].isEmpty();
            assertNoneAgrees(strs, r, none);
            if (!none)
                Assertions.assertEquals(STRS[r], strs.get(r));
        }

        BoolColumn bools = (BoolColumn) df.getColumn("bools");
        Assertions.assertEquals(Types.BOOL, bools.getType());
        for (int r = 0; r < BOOLS.length; r++) {
            assertNoneAgrees(bools, r, false); // bool is never none
            boolean expected = BOOLS[r] != null && BOOLS[r];
            Assertions.assertEquals(expected, bools.get(r));
        }

        FloatColumn floats = (FloatColumn) df.getColumn("floats");
        Assertions.assertEquals(Types.FLOAT, floats.getType());
        Assertions.assertFalse(floats.isDoublePrecision());
        for (int r = 0; r < FLOATS.length; r++) {
            boolean none = FLOATS[r] == null;
            assertNoneAgrees(floats, r, none);
            if (none)
                continue;
            if (Float.isNaN(FLOATS[r]))
                Assertions.assertTrue(Float.isNaN(floats.get(r)));
            else
                Assertions.assertEquals(FLOATS[r], floats.get(r), 0.0F);
        }

        FloatColumn doubles = (FloatColumn) df.getColumn("doubles");
        Assertions.assertEquals(Types.FLOAT, doubles.getType());
        Assertions.assertTrue(doubles.isDoublePrecision());
        for (int r = 0; r < DOUBLES.length; r++) {
            boolean none = DOUBLES[r] == null;
            assertNoneAgrees(doubles, r, none);
            if (none)
                continue;
            if (Double.isNaN(DOUBLES[r]))
                Assertions.assertTrue(Double.isNaN(doubles.getDouble(r)));
            else
                Assertions.assertEquals(DOUBLES[r], doubles.getDouble(r), 0.0);
        }

        DateTimeColumn dates = (DateTimeColumn) df.getColumn("dates");
        Assertions.assertEquals(Types.DATE_TIME, dates.getType());
        for (int r = 0; r < DATES.length; r++) {
            boolean none = DATES[r] == null;
            assertNoneAgrees(dates, r, none);
            if (!none)
                Assertions.assertEquals(DATES[r], dates.get(r), 0.0);
        }

        BigIntColumn bigs = (BigIntColumn) df.getColumn("bigs");
        Assertions.assertEquals(Types.BIG_INT, bigs.getType());
        for (int r = 0; r < BIGS.length; r++) {
            boolean none = BIGS[r] == null || BIGS[r].isEmpty();
            assertNoneAgrees(bigs, r, none);
            if (!none)
                Assertions.assertEquals(BIGS[r], bigs.get(r));
        }
    }

    @Test
    public void testFullFrameRoundTrip() {
        DataFrame df = buildFullFrame();
        byte[] bytes = df.toByteArray();
        DataFrame decoded = DataFrame.fromByteArray(bytes);
        assertFullFrame(decoded);
    }

    @Test
    public void testZeroRowFrame() {
        DataFrame df = new DataFrame();
        df.name = "empty";
        df.addColumn(new IntColumn("i", 0));
        df.addColumn(new StringColumn("s", 0));
        df.addColumn(new DateTimeColumn("d", 0));
        df.addColumn(FloatColumn.double64("f64", new Double[0]));

        DataFrame decoded = DataFrame.fromByteArray(df.toByteArray());
        Assertions.assertEquals("empty", decoded.name);
        Assertions.assertEquals(Integer.valueOf(0), decoded.rowCount);
        Assertions.assertEquals(4, decoded.getColumnCount());
        Assertions.assertEquals(Types.INT, decoded.getColumn("i").getType());
        Assertions.assertEquals(Types.STRING, decoded.getColumn("s").getType());
        Assertions.assertEquals(Types.DATE_TIME, decoded.getColumn("d").getType());
        Assertions.assertEquals(Types.FLOAT, decoded.getColumn("f64").getType());
    }

    @Test
    public void testMultiTableBlob() {
        DataFrame a = new DataFrame();
        a.name = "a";
        a.addColumn(new IntColumn("ints", INTS));

        DataFrame b = new DataFrame();
        b.name = "b";
        b.addColumn(new StringColumn("strs", STRS));

        byte[] bytes = new TablesBlob(new DataFrame[]{a, b}).toByteArray();
        TablesBlob blob = TablesBlob.fromByteArray(bytes);

        DataFrame da = blob.getTable(0);
        DataFrame db = blob.getTable(1);
        Assertions.assertEquals("a", da.name);
        Assertions.assertEquals("b", db.name);

        IntColumn ints = (IntColumn) da.getColumn("ints");
        for (int r = 0; r < INTS.length; r++) {
            boolean none = INTS[r] == null;
            Assertions.assertEquals(none, ints.isNone(r));
            if (!none)
                Assertions.assertEquals(INTS[r], ints.get(r));
        }
        StringColumn strs = (StringColumn) db.getColumn("strs");
        Assertions.assertEquals(EMOJI, strs.get(4));
    }

    @Test
    public void testColumnNameUniquificationSurvives() {
        DataFrame df = new DataFrame();
        df.name = "dup";
        df.addColumn(new IntColumn("x", new Integer[]{1, 2, 3}));
        df.addColumn(new IntColumn("x", new Integer[]{4, 5, 6}));
        // Writer uniquified the second column name.
        Assertions.assertEquals("x", df.getColumn(0).getName());
        Assertions.assertEquals("x (1)", df.getColumn(1).getName());

        DataFrame decoded = DataFrame.fromByteArray(df.toByteArray());
        Assertions.assertEquals("x", decoded.getColumn(0).getName());
        Assertions.assertEquals("x (1)", decoded.getColumn(1).getName());
        Assertions.assertEquals(Integer.valueOf(6), ((IntColumn) decoded.getColumn(1)).get(2));
    }

    // The Java writer only emits int:raw (id 1). These focused tests exercise the
    // int:pattern (id 2) and int:bitIntList (id 4) decoder ports directly, by
    // hand-building their wire form per the Dart serialize layout. Full cross-
    // language coverage (incl. int:rle id 3) lands with the WO-2 Dart fixtures.

    @Test
    public void testIntPatternDecode() {
        // Arithmetic sequence [4,5,6,7,8]: start=4, step=1, blockLength=1,
        // blocksPerCycle=length=5 (IntSequencePattern.serialize order).
        BufferAccessor buf = new BufferAccessor();
        buf.writeInt32(2); // int:pattern encoder id
        buf.writeInt32(4); // start
        buf.writeInt32(1); // step
        buf.writeInt32(1); // blockLength
        buf.writeInt32(5); // blocksPerCycle
        buf.writeInt32(5); // length

        IntColumn col = new IntColumn("p", 0);
        col.decode(new BufferAccessor(buf.toUint8List()));
        Assertions.assertEquals(5, col.getLength());
        for (int i = 0; i < 5; i++)
            Assertions.assertEquals(Integer.valueOf(4 + i), col.get(i));
    }

    @Test
    public void testIntBitIntListDecode() {
        // Values [0,1,2,3,None] with min=0. Element code = None ? 0 : v - min + 1,
        // packed MSB-first, 3 bits each (mirrors BitIntList._setInt / serialize).
        int bits = 3;
        int min = 0;
        int[] values = {0, 1, 2, 3, IntColumn.None};
        int intsPer32 = 32 / bits;
        int[] words = new int[(values.length + intsPer32 - 1) / intsPer32];
        for (int i = 0; i < values.length; i++) {
            int code = values[i] == IntColumn.None ? 0 : values[i] - min + 1;
            setBits(words, i, bits, code);
        }

        BufferAccessor buf = new BufferAccessor();
        buf.writeInt32(4);          // int:bitIntList encoder id
        buf.writeInt32(bits);       // bits
        buf.writeInt64(min);        // min
        buf.writeInt64(values.length); // length
        buf.writeInt8((byte)0);     // archive
        buf.writeUint32List(words);

        IntColumn col = new IntColumn("b", 0);
        col.decode(new BufferAccessor(buf.toUint8List()));
        Assertions.assertEquals(values.length, col.getLength());
        for (int i = 0; i < values.length; i++)
            Assertions.assertEquals(Integer.valueOf(values[i]), col.get(i));
        Assertions.assertTrue(col.isNone(4));
    }

    // Reference bit packer mirroring Dart BitIntList._setInt (MSB-first per word).
    private static void setBits(int[] data, int i, int bits, int value) {
        long mask = (1L << bits) - 1;
        int intsPer32 = 32 / bits;
        long ui32 = data[i / intsPer32] & 0xFFFFFFFFL;
        int shift = bits * (intsPer32 - i % intsPer32 - 1);
        ui32 = ui32 & ~(mask << shift);
        ui32 = ui32 | ((value & mask) << shift);
        data[i / intsPer32] = (int) ui32;
    }

    @Test
    public void testTruncatedBufferThrowsCleanException() {
        byte[] bytes = buildFullFrame().toByteArray();
        byte[] truncated = new byte[bytes.length / 2];
        System.arraycopy(bytes, 0, truncated, 0, truncated.length);
        RuntimeException e = Assertions.assertThrows(RuntimeException.class,
                () -> DataFrame.fromByteArray(truncated));
        Assertions.assertFalse(e instanceof IndexOutOfBoundsException, "must not leak a raw OOB");
    }

    // A truncated multi-table blob must fail cleanly from TablesBlob.fromByteArray +
    // getTable(0), not leak an ArrayIndexOutOfBoundsException with a garbage index.
    @Test
    public void testTruncatedMultiTableBlobThrowsCleanException() {
        DataFrame a = new DataFrame();
        a.name = "a";
        a.addColumn(new IntColumn("ints", INTS));
        DataFrame b = new DataFrame();
        b.name = "b";
        b.addColumn(new StringColumn("strs", STRS));

        byte[] bytes = new TablesBlob(new DataFrame[]{a, b}).toByteArray();
        byte[] truncated = new byte[bytes.length / 2];
        System.arraycopy(bytes, 0, truncated, 0, truncated.length);

        RuntimeException e = Assertions.assertThrows(RuntimeException.class,
                () -> TablesBlob.fromByteArray(truncated).getTable(0));
        Assertions.assertFalse(e instanceof IndexOutOfBoundsException, "must not leak a raw OOB");
    }

    // Guards the highest-risk change in WO-1: a default float column must encode as
    // id 1 (raw32) and round-trip single-precision; only double64() emits id 5.
    // A regression flipping the default to id 5 would silently alter query-result
    // serialization platform-wide.
    @Test
    public void testFloatDefaultEncoderPathIsId1() {
        FloatColumn f = new FloatColumn("f", new Float[]{1.5F, -2.25F, null});
        BufferAccessor buf = new BufferAccessor();
        f.encode(buf);
        Assertions.assertEquals(1, new BufferAccessor(buf.toUint8List()).readInt32(),
                "default float column must encode as id 1 (raw32)");

        FloatColumn d = FloatColumn.double64("d", new Double[]{1.5, -2.25});
        BufferAccessor buf2 = new BufferAccessor();
        d.encode(buf2);
        Assertions.assertEquals(5, new BufferAccessor(buf2.toUint8List()).readInt32(),
                "double64 column must encode as id 5 (raw64)");

        // Default float round-trips single precision.
        DataFrame df = new DataFrame();
        df.name = "f";
        df.addColumn(new FloatColumn("f", new Float[]{1.5F, -2.25F, null}));
        FloatColumn decoded = (FloatColumn) DataFrame.fromByteArray(df.toByteArray()).getColumn("f");
        Assertions.assertFalse(decoded.isDoublePrecision());
        Assertions.assertEquals(1.5F, decoded.get(0), 0.0F);
        Assertions.assertTrue(decoded.isNone(2));
    }

    // readUint16List must be unsigned: values >= 0x8000 read back as 32768..65535,
    // not sign-extended negatives.
    @Test
    public void testReadUint16ListIsUnsigned() {
        BufferAccessor buf = new BufferAccessor();
        buf.writeInt16((short) BufferAccessor.TYPE_UINT_16_LIST); // type code
        buf.writeInt64(3);                                        // count
        buf.writeInt16((short) 0xFFFF);
        buf.writeInt16((short) 0x8000);
        buf.writeInt16((short) 0x1234);

        int[] r = new BufferAccessor(buf.toUint8List()).readUint16List();
        Assertions.assertArrayEquals(new int[]{0xFFFF, 0x8000, 0x1234}, r);
    }
}
