package serialization;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Random;

// IntColumn.encode picks the encoder exactly like Dart's IntMeta / getBestEncodingEstimate.
public class IntColumnEncoderSelectionTest {
    private static final int N = IntColumn.None;

    @AfterEach
    public void restoreSwitch() {
        IntColumn.ADVANCED_ENCODERS = true;
    }

    static byte[] encode(Column<?> col) {
        BufferAccessor buf = new BufferAccessor();
        col.encode(buf);
        return buf.toUint8List();
    }

    static int encoderId(byte[] bytes) {
        return new BufferAccessor(bytes).readInt32();
    }

    // Encoder id of the nested index column of a string payload (id 0 or 4).
    static int indexEncoderId(StringColumn col) {
        BufferAccessor buf = new BufferAccessor(encode(col));
        int id = buf.readInt32();
        if (id == 4)
            serialization.codecs.StringTokens.decode(buf, "");
        else {
            Assertions.assertEquals(0, id);
            buf.readStringList();
        }
        return buf.readInt32();
    }

    private static IntColumn decode(byte[] bytes) {
        IntColumn col = new IntColumn("", 0);
        col.decode(new BufferAccessor(bytes));
        return col;
    }

    private static void assertEncodes(int[] values, int expectedId) {
        IntColumn col = new IntColumn("c", values);
        byte[] bytes = encode(col);
        Assertions.assertEquals(expectedId, encoderId(bytes), "encoder id");
        Assertions.assertArrayEquals(values, (int[]) decode(bytes).toArray());
    }

    static int[] sequential(int count) {
        int[] values = new int[count];
        for (int i = 0; i < count; i++)
            values[i] = i + 1;
        return values;
    }

    // Runs of irregular length and occasional skipped values: zero deltas with rare
    // jumps (rle territory), while a regular i / groupSize would be an int:pattern.
    static int[] grouped(int count, int groupSize) {
        int[] values = new int[count];
        int g = 0;
        int left = groupSize;
        for (int i = 0; i < count; i++) {
            if (left == 0) {
                g += g % 5 == 4 ? 2 : 1;
                left = groupSize + g % 11;
            }
            values[i] = g;
            left--;
        }
        return values;
    }

    static String[] plates(int[] groups) {
        String[] plates = new String[groups.length];
        for (int i = 0; i < groups.length; i++)
            plates[i] = String.format("PLT-%04d", groups[i]);
        return plates;
    }

    static int[] smallRange(int count, int bound, long seed) {
        Random rnd = new Random(seed);
        int[] values = new int[count];
        for (int i = 0; i < count; i++)
            values[i] = i % 37 == 0 ? N : rnd.nextInt(bound);
        return values;
    }

    static int[] highEntropy(int count, long seed) {
        Random rnd = new Random(seed);
        int[] values = new int[count];
        for (int i = 0; i < count; i++) {
            values[i] = rnd.nextInt();
            if (values[i] == N)
                values[i] = 0;
        }
        return values;
    }

    @Test
    public void sequentialPicksPatternWithZeroPayload() {
        byte[] bytes = encode(new IntColumn("seq", sequential(100000)));
        Assertions.assertEquals(2, encoderId(bytes));
        Assertions.assertEquals(24, bytes.length); // id + 5 pattern ints, no archive byte
        Assertions.assertArrayEquals(sequential(100000), (int[]) decode(bytes).toArray());
    }

    @Test
    public void groupedPicksRle() {
        assertEncodes(grouped(100000, 384), 3);
        assertEncodes(grouped(1000, 50), 3);
    }

    @Test
    public void smallRangePicksBitIntList() {
        // 4-bit values pack identically in bitIntList and bitPacked; the tie keeps the legacy id.
        assertEncodes(smallRange(100000, 10, 1), 4);
        // 5-bit values waste 2 bits/word in bitIntList; dense int:bitPacked (GROK-20761) wins.
        assertEncodes(smallRange(1000, 24, 2), 5);
    }

    @Test
    public void highEntropyStaysRaw() {
        assertEncodes(highEntropy(100000, 3), 1);
    }

    // Dart: a constant None run of 3+ is an int:pattern (step 0); shorter all-None and empty columns stay raw.
    @Test
    public void emptyOrAllNoneFollowsDart() {
        assertEncodes(new int[0], 1);
        assertEncodes(new int[]{N}, 1);
        assertEncodes(new int[]{N, N}, 1);
        assertEncodes(new int[]{N, N, N, N}, 2);
    }

    // Constants shorter than a pattern: bitIntList is 1 bit per element (4 bytes per 32), so a single
    // value ties raw and stays raw (Dart's msb(max - min) = 0 estimate would pick bitIntList there).
    @Test
    public void regularBlocksPickPatternAndConstantsPickBitIntList() {
        assertEncodes(new int[]{42}, 1);
        assertEncodes(new int[]{7, 7}, 4);
        assertEncodes(new int[]{7, 7, 7}, 2);
        int[] regular = new int[100000];
        for (int i = 0; i < regular.length; i++)
            regular[i] = i / 384;
        assertEncodes(regular, 2);
    }

    @Test
    public void switchOffRestoresRaw() {
        IntColumn.ADVANCED_ENCODERS = false;
        assertEncodes(sequential(1000), 1);
        assertEncodes(grouped(1000, 50), 1);
        assertEncodes(smallRange(1000, 10, 1), 1);
    }

    @Test
    public void stringIndexColumnUsesBestEncoder() {
        int[] regular = new int[100000];
        for (int i = 0; i < regular.length; i++)
            regular[i] = i / 384;
        Assertions.assertEquals(2, indexEncoderId(new StringColumn("plate", plates(regular))));
        Assertions.assertEquals(3, indexEncoderId(new StringColumn("plate", plates(grouped(100000, 384)))));

        String[] few = new String[1000];
        String[] cats = {"alpha", "beta", "gamma"};
        Random rnd = new Random(5);
        for (int i = 0; i < few.length; i++)
            few[i] = cats[rnd.nextInt(3)];
        Assertions.assertEquals(4, indexEncoderId(new StringColumn("few", few)));

        String[] same = new String[1000];
        for (int i = 0; i < same.length; i++)
            same[i] = "only";
        Assertions.assertEquals(2, indexEncoderId(new StringColumn("same", same)));
    }

    @Test
    public void bigIntIntModeUsesBestEncoder() {
        BigIntColumn col = new BigIntColumn("ids", 10);
        col.setDowncastAllowed(true);
        for (int i = 1; i <= 100000; i++)
            col.add((long) i);
        byte[] bytes = encode(col);
        Assertions.assertEquals(Types.INT, col.getType());
        Assertions.assertEquals(2, encoderId(bytes));
        Assertions.assertArrayEquals(sequential(100000), (int[]) decode(bytes).toArray());
    }
}
