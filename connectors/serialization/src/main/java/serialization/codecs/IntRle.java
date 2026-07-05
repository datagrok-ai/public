package serialization.codecs;

import serialization.BufferAccessor;
import serialization.IntColumn;

// Decode-only port of ddt/lib/src/serialization/int_rle.dart:214-248.
// Run-length encoding with delta encoding pre-processing.
public class IntRle {
    // Block flag mask.
    static final int BLOCK_FLAG_MASK = 0x8000;

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
