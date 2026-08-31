package serialization.codecs;

import serialization.BufferAccessor;
import serialization.ByteData;
import serialization.ColumnEncoderArchiveType;
import serialization.IntColumn;
import serialization.Zlib;

// Port of ddt/lib/src/serialization/bit_int_list.dart.
// An efficient storage of ints where each element takes the number of bits
// necessary to represent (max - min). Element mapping: raw==0 -> None,
// raw v -> v + min - 1. Bit packing is MSB-first within each 32-bit word.
public class BitIntList {
    public static final int None = IntColumn.None;

    private int[] data;
    private int bits = 32;
    private long min = 0;
    private int length = 0;

    private BitIntList() {
    }

    public static BitIntList fromList(int[] values, int start, int count) {
        long min = 0;
        long max = 0;
        boolean any = false;
        for (int i = start; i < start + count; i++) {
            if (values[i] == None)
                continue;
            min = any ? Math.min(min, values[i]) : values[i];
            max = any ? Math.max(max, values[i]) : values[i];
            any = true;
        }

        BitIntList b = new BitIntList();
        b.length = count;
        b.min = min;
        b.bits = count == 0 ? 0 : msb(max - min + 1);
        b.data = new int[b.bits == 0 ? 0 : wordCount(count, b.bits)];
        for (int i = 0; i < count; i++)
            b.set(i, values[start + i]);
        return b;
    }

    // What fromList writes; bit_int_list.dart:52 matches since GROK-20761.
    public static int sizeInBytes(int count, long min, long max) {
        int bitsPerInt = count == 0 ? 0 : msb(max - min + 1);
        return (bitsPerInt == 0 ? 0 : wordCount(count, bitsPerInt)) * 4;
    }

    public static int msb(long v) {
        return v <= 0 ? 0 : 64 - Long.numberOfLeadingZeros(v);
    }

    private static int wordCount(int count, int bits) {
        int intsPer32 = 32 / bits;
        return (count + intsPer32 - 1) / intsPer32;
    }

    public static BitIntList fromBuffer(BufferAccessor buffer) {
        BitIntList b = new BitIntList();
        b.bits = buffer.readInt32();
        b.min = buffer.readInt64();
        b.length = (int) buffer.readInt64();

        if (buffer.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB) {
            byte[] bytes = Zlib.inflate(buffer.readUint8List());
            b.data = ByteData.toUint32List(bytes);
        } else
            b.data = buffer.readUint32List();

        return b;
    }

    public void serialize(BufferAccessor buffer) {
        buffer.writeInt32(bits);
        buffer.writeInt64(min);
        buffer.writeInt64(length);
        buffer.writeInt8((byte) ColumnEncoderArchiveType.ARCHIVE_TYPE_NONE);
        buffer.writeUint32List(data);
    }

    public int get(int idx) {
        if (idx < 0 || idx >= length)
            throw new RuntimeException("Index out of bounds: " + idx);

        if (bits == 0)
            return (int) min;

        long val = getInt(data, idx, bits);
        return val == 0 ? None : (int) (val + min - 1);
    }

    private void set(int idx, int value) {
        if (idx < 0 || idx >= length)
            throw new RuntimeException("Index out of bounds: " + idx);

        if (bits != 0)
            setInt(data, idx, bits, value == None ? 0 : value - min + 1);
    }

    public int[] toInt32List() {
        int[] r = new int[length];
        for (int i = 0; i < length; i++)
            r[i] = get(i);
        return r;
    }

    private static long getInt(int[] data, int i, int bits) {
        long mask = (1L << bits) - 1;
        int intsPer32 = 32 / bits;
        long ui32 = data[i / intsPer32] & 0xFFFFFFFFL;
        long x = ui32 >>> (bits * (intsPer32 - i % intsPer32 - 1));
        return x & mask;
    }

    private static void setInt(int[] data, int i, int bits, long value) {
        long mask = (1L << bits) - 1;
        int intsPer32 = 32 / bits;
        long ui32 = data[i / intsPer32] & 0xFFFFFFFFL;
        int shift = bits * (intsPer32 - i % intsPer32 - 1);
        ui32 = ui32 & ~(mask << shift);
        ui32 = ui32 | (value << shift);
        data[i / intsPer32] = (int) ui32;
    }
}
