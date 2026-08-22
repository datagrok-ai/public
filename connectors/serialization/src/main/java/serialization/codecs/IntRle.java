package serialization.codecs;

import serialization.BufferAccessor;
import serialization.IntColumn;

// Port of ddt/lib/src/serialization/int_rle.dart.
// Run-length encoding with delta encoding pre-processing.
public class IntRle {
    static final int MAX_BLOCK_LENGTH = 32767;
    static final int BLOCK_FLAG_MASK = 0x8000;

    private int lengthRepeats;
    private int lengthSymbols;
    private boolean estimated;

    // Deltas are 64-bit like Dart ints; symbols truncate to int32 as Dart's Int32List does
    // and decode's int32 accumulation restores the values mod 2^32. Returns -1 when not encodable.
    public int estimate(int[] values, int length, int minimum) {
        if (length < 3)
            return -1;

        int none = minimum - 1;
        long vn = replaceNone(values[1], none);
        long vnm1 = replaceNone(values[0], none);
        long prev = vn - vnm1;
        long curr;
        int cntSymbols = 0;
        boolean isBlock = false;
        long max = prev;
        long min = prev;
        lengthRepeats = 0;
        lengthSymbols = 0;

        for (int n = 1; n < length; n++) {
            vn = replaceNone(values[n], none);
            curr = vn - vnm1;
            vnm1 = vn;

            if (curr < min)
                min = curr;

            if (curr > max)
                max = curr;

            if (curr == prev) {
                if (!isBlock) {
                    cntSymbols++;
                    if (cntSymbols > MAX_BLOCK_LENGTH) {
                        lengthRepeats++;
                        lengthSymbols++;
                        cntSymbols = 1;
                    }
                } else {
                    if (cntSymbols > 0)
                        lengthRepeats++;
                    cntSymbols = 2;
                    isBlock = false;
                }
            } else {
                if (isBlock) {
                    lengthSymbols++;
                    cntSymbols++;
                    if (cntSymbols > MAX_BLOCK_LENGTH) {
                        lengthRepeats++;
                        cntSymbols = 1;
                    }
                } else {
                    lengthRepeats++;
                    lengthSymbols++;
                    cntSymbols = 0;
                    isBlock = true;
                }
            }

            prev = curr;
        }

        if (!isBlock) {
            lengthRepeats++;
            lengthSymbols++;
        } else {
            cntSymbols++;
            if (cntSymbols > MAX_BLOCK_LENGTH)
                lengthRepeats++;
            lengthRepeats++;
            lengthSymbols++;
        }

        if (BitIntList.msb(max - min) > 31)
            return -1;

        estimated = true;
        // Repeat entries are 1..0xFFFF, so 16 bits each on the wire (Dart's estimate shape, with sizeInBytes's +1).
        int repeatsSizeInBytes = BitIntList.sizeInBytes(lengthRepeats, 1, 0xFFFF);
        int symbolsSizeInBytes = BitIntList.sizeInBytes(lengthSymbols, min, max);
        return repeatsSizeInBytes + symbolsSizeInBytes;
    }

    public void encode(BufferAccessor buffer, int[] values, int length, int minimum) {
        if (!estimated)
            throw new IllegalStateException("Perform 'estimate' first.");

        int[] repeats = new int[lengthRepeats];
        int[] symbols = new int[lengthSymbols];

        int none = minimum - 1;
        long vn = replaceNone(values[1], none);
        long vnm1 = replaceNone(values[0], none);
        long prev = vn - vnm1;
        long curr;
        int idxSymbols = 0;
        int idxRepeats = 0;
        int cntSymbols = 0;
        boolean isBlock = false;

        for (int n = 1; n < length; n++) {
            vn = replaceNone(values[n], none);
            curr = vn - vnm1;
            vnm1 = vn;

            if (curr == prev) {
                if (!isBlock) {
                    cntSymbols++;
                    if (cntSymbols > MAX_BLOCK_LENGTH) {
                        repeats[idxRepeats++] = MAX_BLOCK_LENGTH;
                        symbols[idxSymbols++] = (int) curr;
                        cntSymbols = 1;
                    }
                } else {
                    if (cntSymbols > 0)
                        repeats[idxRepeats++] = cntSymbols | BLOCK_FLAG_MASK;
                    cntSymbols = 2;
                    isBlock = false;
                }
            } else {
                if (isBlock) {
                    symbols[idxSymbols++] = (int) prev;
                    cntSymbols++;
                    if (cntSymbols > MAX_BLOCK_LENGTH) {
                        repeats[idxRepeats++] = MAX_BLOCK_LENGTH | BLOCK_FLAG_MASK;
                        cntSymbols = 1;
                    }
                } else {
                    repeats[idxRepeats++] = cntSymbols;
                    symbols[idxSymbols++] = (int) prev;
                    cntSymbols = 0;
                    isBlock = true;
                }
            }

            prev = curr;
        }

        if (!isBlock) {
            repeats[idxRepeats] = cntSymbols;
            symbols[idxSymbols] = (int) prev;
        } else {
            cntSymbols++;
            if (cntSymbols > MAX_BLOCK_LENGTH)
                repeats[idxRepeats++] = MAX_BLOCK_LENGTH | BLOCK_FLAG_MASK;
            repeats[idxRepeats] = cntSymbols | BLOCK_FLAG_MASK;
            symbols[idxSymbols] = (int) prev;
        }

        buffer.writeInt32(replaceNone(values[0], none));
        buffer.writeInt32(none);
        buffer.writeInt32(length);
        BitIntList.fromList(symbols, 0, symbols.length).serialize(buffer);
        BitIntList.fromList(repeats, 0, repeats.length).serialize(buffer);

        estimated = false;
    }

    private static int replaceNone(int x, int none) {
        return x == IntColumn.None ? none : x;
    }

    public static int[] decode(BufferAccessor buffer) {
        int prev = buffer.readInt32();
        int none = buffer.readInt32();
        int[] decoded = new int[buffer.readInt32()];
        int[] symbols = BitIntList.fromBuffer(buffer).toInt32List();
        int[] repeats = BitIntList.fromBuffer(buffer).toInt32List();

        int idxDecoded = 1;
        int acc = prev;
        int idxSymbols = 0;
        decoded[0] = prev;
        for (int n = 0; n < repeats.length; n++) {
            int cnt = repeats[n] & (~BLOCK_FLAG_MASK);
            boolean isBlock = (repeats[n] & BLOCK_FLAG_MASK) == BLOCK_FLAG_MASK;
            if (isBlock) {
                while (cnt-- > 0) {
                    acc += symbols[idxSymbols++];
                    decoded[idxDecoded++] = acc;
                }
            } else {
                while (cnt-- > 0) {
                    acc += symbols[idxSymbols];
                    decoded[idxDecoded++] = acc;
                }
                idxSymbols++;
            }
        }

        // Restore Nones
        for (int n = 0; n < decoded.length; n++)
            if (decoded[n] == none)
                decoded[n] = IntColumn.None;

        return decoded;
    }
}
