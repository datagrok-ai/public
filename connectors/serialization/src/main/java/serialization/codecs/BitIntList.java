package serialization.codecs;

import serialization.BufferAccessor;
import serialization.ByteData;
import serialization.ColumnEncoderArchiveType;
import serialization.IntColumn;
import serialization.Zlib;

// Decode-only port of ddt/lib/src/serialization/bit_int_list.dart:57-147.
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

    public int get(int idx) {
        if (idx < 0 || idx >= length)
            throw new RuntimeException("Index out of bounds: " + idx);

        if (bits == 0)
            return (int) min;

        long val = getInt(data, idx, bits);
        return val == 0 ? None : (int) (val + min - 1);
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
}
