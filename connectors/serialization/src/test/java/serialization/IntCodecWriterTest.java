package serialization;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import serialization.codecs.BitIntList;
import serialization.codecs.IntRle;
import serialization.codecs.IntSequencePattern;

import java.util.Random;

// Write side of the int codecs (ports of bit_int_list.dart, int_pattern.dart,
// int_rle.dart) round-tripped through the pre-existing Java decoders.
public class IntCodecWriterTest {
    private static final int N = IntColumn.None;

    private static int[] randomWithNones(Random rnd, int count, int min, int max, double noneShare) {
        int[] values = new int[count];
        for (int i = 0; i < count; i++)
            values[i] = rnd.nextDouble() < noneShare ? N : (int) (min + (long) (rnd.nextDouble() * ((long) max - min + 1)));
        return values;
    }

    private static byte[] bitIntListBytes(int[] values) {
        BufferAccessor buf = new BufferAccessor();
        BitIntList.fromList(values, 0, values.length).serialize(buf);
        return buf.toUint8List();
    }

    private static void assertBitIntListRoundTrip(int[] values) {
        int[] decoded = BitIntList.fromBuffer(new BufferAccessor(bitIntListBytes(values))).toInt32List();
        Assertions.assertArrayEquals(values, decoded);
    }

    @Test
    public void bitIntListRoundTripsRandomArraysWithNones() {
        Random rnd = new Random(42);
        assertBitIntListRoundTrip(randomWithNones(rnd, 1000, 0, 9, 0.1));
        assertBitIntListRoundTrip(randomWithNones(rnd, 1000, -5, 5, 0.1));
        assertBitIntListRoundTrip(randomWithNones(rnd, 777, 100000, 100127, 0.0));
        assertBitIntListRoundTrip(randomWithNones(rnd, 333, -70000, 70000, 0.3));
        assertBitIntListRoundTrip(randomWithNones(rnd, 100, Integer.MIN_VALUE + 1, Integer.MAX_VALUE, 0.1));
        assertBitIntListRoundTrip(new int[]{7, 7, 7});
        assertBitIntListRoundTrip(new int[]{N, N, N});
        assertBitIntListRoundTrip(new int[]{N, 3});
        assertBitIntListRoundTrip(new int[]{42});
        assertBitIntListRoundTrip(new int[0]);
    }

    @Test
    public void bitIntListPacksMsbFirstPerWord() {
        // Values [0,1,2,3,None] with min=0 -> codes 1,2,3,4,0 in 3 bits, MSB-first:
        // 001 010 011 100 000 ... -> 0x0A700000 (matches D42ReaderRoundTripTest.setBits).
        byte[] bytes = bitIntListBytes(new int[]{0, 1, 2, 3, N});
        BufferAccessor buf = new BufferAccessor(bytes);
        Assertions.assertEquals(3, buf.readInt32());
        Assertions.assertEquals(0, buf.readInt64());
        Assertions.assertEquals(5, buf.readInt64());
        Assertions.assertEquals(0, buf.readInt8());
        Assertions.assertArrayEquals(new int[]{0x0A700000}, buf.readUint32List());
    }

    // sizeInBytes is what fromList writes: msb(max - min + 1) bits per element (Dart estimates msb(max - min)).
    @Test
    public void sizeInBytesMatchesFromList() {
        Assertions.assertEquals(0, BitIntList.sizeInBytes(0, 0, 100));
        Assertions.assertEquals(4, BitIntList.sizeInBytes(10, 5, 5));
        Assertions.assertEquals(8, BitIntList.sizeInBytes(32, 0, 1));
        Assertions.assertEquals(12, BitIntList.sizeInBytes(33, 0, 1));
        Assertions.assertEquals(4 * 125, BitIntList.sizeInBytes(1000, 0, 9));
        Assertions.assertEquals(4 * 1000, BitIntList.sizeInBytes(1000, Integer.MIN_VALUE + 1, Integer.MAX_VALUE));
        // 21-byte bitIntList header + 10-byte uint32 list header precede the packed words.
        int[] values = randomWithNones(new Random(1), 1000, 0, 9, 0.0);
        Assertions.assertEquals(BitIntList.sizeInBytes(1000, 0, 9), bitIntListBytes(values).length - 21 - 10);
        // Range 0..7 is a power-of-two range: msb(8) = 4 bits on the wire, where Dart's msb(7) = 3 would undercount.
        values = randomWithNones(new Random(2), 1000, 0, 7, 0.0);
        Assertions.assertEquals(4 * 125, BitIntList.sizeInBytes(1000, 0, 7));
        Assertions.assertEquals(4 * 125, bitIntListBytes(values).length - 21 - 10);
        values = new int[]{5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
        Assertions.assertEquals(BitIntList.sizeInBytes(10, 5, 5), bitIntListBytes(values).length - 21 - 10);
    }

    private static void assertPattern(int[] values) {
        IntSequencePattern pattern = IntSequencePattern.fromList(values, values.length);
        Assertions.assertNotNull(pattern, "pattern expected");
        BufferAccessor buf = new BufferAccessor();
        pattern.serialize(buf);
        Assertions.assertEquals(20, buf.toUint8List().length);
        int[] decoded = new IntSequencePattern(new BufferAccessor(buf.toUint8List())).toInt32List();
        Assertions.assertArrayEquals(values, decoded);
    }

    @Test
    public void patternRoundTrips() {
        assertPattern(new int[]{4, 5, 6, 7, 8});
        assertPattern(new int[]{5, 7, 9, 11});
        assertPattern(new int[]{1, 1, 2, 2, 3, 3, 4, 4});
        assertPattern(new int[]{10, 10, 10, 10});
        assertPattern(new int[]{N, N, N, N});
        assertPattern(new int[]{1, 2, 3, 1, 2, 3, 1});
        assertPattern(new int[]{5, 5, 6, 6, 5, 5, 6, 6, 5});
        assertPattern(new int[]{100, 90, 80, 70});
        assertPattern(new int[]{1, 1, 2, 2, 3}); // a truncated last block still matches
        int[] seq = new int[100000];
        for (int i = 0; i < seq.length; i++)
            seq[i] = i + 1;
        assertPattern(seq);
        int[] cycle = new int[1000];
        for (int i = 0; i < cycle.length; i++)
            cycle[i] = 1 + (i % 24) / 3;
        assertPattern(cycle);
        // A step exceeding int32 wraps on the wire but decodes mod 2^32 on both sides.
        assertPattern(new int[]{-2000000000, 2000000000, -2000000000, 2000000000});
        // No cycle: blocksPerCycle == length, so blocksPerCycle * blockLength exceeds int32.
        assertPattern(blocks(2, 65536));
        assertPattern(blocks(3, 65536));
        assertPattern(blocks(2, 46341));
    }

    private static int[] blocks(int count, int blockLength) {
        int[] values = new int[count * blockLength];
        for (int i = 0; i < values.length; i++)
            values[i] = i / blockLength;
        return values;
    }

    @Test
    public void patternRejectsNonPatterns() {
        Assertions.assertNull(IntSequencePattern.fromList(new int[]{1, 2}, 2));
        Assertions.assertNull(IntSequencePattern.fromList(new int[]{1, 2, 4}, 3));
        Assertions.assertNull(IntSequencePattern.fromList(new int[]{1, 1, 2, 2, 3, 4}, 6));
        Assertions.assertNull(IntSequencePattern.fromList(new int[]{1, 2, 3, 1, 2, 4}, 6));
        Assertions.assertNull(IntSequencePattern.fromList(new int[]{3, 3, 3, 4, 5}, 5));
        // Wrapped int32 delta must not be mistaken for a real one: 2e9 - (-2e9) != -294967296.
        Assertions.assertNull(IntSequencePattern.fromList(new int[]{-2000000000, 2000000000, 1705032704}, 3));
    }

    private static int min(int[] values) {
        int min = Integer.MAX_VALUE;
        for (int v : values)
            if (v != N)
                min = Math.min(min, v);
        return min;
    }

    private static void assertRle(int[] values) {
        IntRle rle = new IntRle();
        int size = rle.estimate(values, values.length, min(values));
        Assertions.assertTrue(size > 0, "rle expected to be applicable");
        BufferAccessor buf = new BufferAccessor();
        rle.encode(buf, values, values.length, min(values));
        int[] decoded = IntRle.decode(new BufferAccessor(buf.toUint8List()));
        Assertions.assertArrayEquals(values, decoded);
    }

    @Test
    public void rleRoundTrips() {
        Random rnd = new Random(7);
        int[] walk = new int[5000];
        int v = 0;
        for (int i = 0; i < walk.length; i++) {
            if (rnd.nextInt(10) == 0)
                walk[i] = N;
            else {
                v += rnd.nextInt(3) - 1;
                walk[i] = v;
            }
        }
        assertRle(walk);

        int[] grouped = new int[10000];
        for (int i = 0; i < grouped.length; i++)
            grouped[i] = i / 50;
        assertRle(grouped);

        // Runs longer than MAX_BLOCK_LENGTH split into several repeat entries.
        int[] constant = new int[70000];
        for (int i = 0; i < constant.length; i++)
            constant[i] = 3;
        assertRle(constant);

        int[] noisy = new int[70000];
        for (int i = 0; i < noisy.length; i++)
            noisy[i] = rnd.nextInt(1000);
        assertRle(noisy);

        // none = min - 1 collides with the None sentinel itself when min is MinValue + 1.
        assertRle(new int[]{Integer.MIN_VALUE + 1, N, Integer.MIN_VALUE + 1, Integer.MIN_VALUE + 2, N});
        assertRle(new int[]{N, N, N, 5, 5, 5});
        assertRle(new int[]{1, 2, 3});
    }

    @Test
    public void rleEstimateRejectsShortOrWideInputs() {
        Assertions.assertEquals(-1, new IntRle().estimate(new int[]{1, 2}, 2, 1));
        Assertions.assertEquals(-1, new IntRle().estimate(new int[]{0, 2000000000, -2000000000, 0}, 4, -2000000000));
        Assertions.assertThrows(IllegalStateException.class,
                () -> new IntRle().encode(new BufferAccessor(), new int[]{1, 2, 3}, 3, 1));
    }
}
