package serialization.codecs;

import serialization.BufferAccessor;
import serialization.IntColumn;

// Port of ddt/lib/src/serialization/bit_packing.dart (GROK-20761); same structure and names.
// Dense frame-of-reference bit packing: value i occupies bits [i*bits, (i+1)*bits) of a
// Uint32List, LSB-first. Codes are (v - base); when the data contains None, code 0 is None
// and other values are (v - base + 1).
// Layout: Int32 bits, Int32 base, Int32 length, Int8 flags (1 = has None), Uint32List words.
public class BitPacking {
    public static final int MAX_BITS = 31;
    public static final int HEADER_BYTES = 23;

    // Returns [min, max, hasNone] of the values, or null when every value is None.
    public static long[] minMax(int[] values, int start, int count) {
        long min = 0;
        long max = 0;
        boolean any = false;
        long hasNone = 0;
        for (int i = start; i < start + count; i++) {
            int v = values[i];
            if (v == IntColumn.None) {
                hasNone = 1;
                continue;
            }
            min = any ? Math.min(min, v) : v;
            max = any ? Math.max(max, v) : v;
            any = true;
        }
        return any ? new long[]{min, max, hasNone} : null;
    }

    // Bits per value, or -1 when the range does not fit MAX_BITS.
    public static int bitsFor(long[] minMax) {
        if (minMax == null)
            return 0;
        int bits = BitIntList.msb(minMax[1] - minMax[0] + minMax[2]);
        return bits > MAX_BITS ? -1 : bits;
    }

    public static int wordCount(int count, int bits) {
        return (int) (((long) count * bits + 31) / 32);
    }

    public static int sizeInBytes(int count, int bits) {
        return HEADER_BYTES + 4 * wordCount(count, bits);
    }

    public static void write(BufferAccessor buf, int[] values, int start, int count) {
        write(buf, values, start, count, minMax(values, start, count));
    }

    public static void write(BufferAccessor buf, int[] values, int start, int count, long[] minMax) {
        int bits = bitsFor(minMax);
        if (bits < 0)
            throw new IllegalArgumentException("BitPacking: value range exceeds " + MAX_BITS + " bits");
        int base = minMax == null ? IntColumn.None : (int) minMax[0];
        int noneShift = minMax == null ? 1 : (int) minMax[2];
        int[] words = new int[wordCount(count, bits)];
        for (int i = 0; bits > 0 && i < count; i++) {
            int v = values[start + i];
            long code = noneShift == 1 && v == IntColumn.None ? 0 : (long) v - base + noneShift;
            long pos = (long) i * bits;
            int w = (int) (pos >> 5);
            int off = (int) (pos & 31);
            words[w] |= (int) (code << off);
            if (off + bits > 32)
                words[w + 1] |= (int) (code >>> (32 - off));
        }
        buf.writeInt32(bits);
        buf.writeInt32(base);
        buf.writeInt32(count);
        buf.writeInt8((byte) noneShift);
        buf.writeUint32List(words);
    }

    public static int[] read(BufferAccessor buf) {
        int bits = buf.readInt32();
        int base = buf.readInt32();
        int count = buf.readInt32();
        int noneShift = buf.readInt8();
        int[] words = buf.readUint32List();
        int[] result = new int[count];
        if (bits == 0) {
            for (int i = 0; i < count; i++)
                result[i] = base;
            return result;
        }
        long mask = (1L << bits) - 1;
        for (int i = 0; i < count; i++) {
            long pos = (long) i * bits;
            int w = (int) (pos >> 5);
            int off = (int) (pos & 31);
            long code = ((words[w] & 0xFFFFFFFFL) >>> off) & mask;
            if (off + bits > 32)
                code |= (words[w + 1] & ((1L << (off + bits - 32)) - 1)) << (32 - off);
            result[i] = noneShift == 1 && code == 0 ? IntColumn.None : (int) (code + base - noneShift);
        }
        return result;
    }
}
