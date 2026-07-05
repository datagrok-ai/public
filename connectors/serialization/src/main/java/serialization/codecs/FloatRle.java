package serialization.codecs;

import serialization.BufferAccessor;
import serialization.ByteData;
import serialization.Zlib;

// Decode-only port of ddt/lib/src/serialization/float_rle.dart:190-259.
// Run-length encoding for float lists. Repeats vector is uint8 (values 0..255);
// bit 7 (BLOCK_FLAG_MASK) marks a symbol block vs a repeated single symbol.
public class FloatRle {
    static final int BLOCK_FLAG_MASK = 0x80;

    public static float[] decode32(BufferAccessor buf) {
        float[] decoded = new float[buf.readInt32()];
        float[] symbols;
        int[] repeats;
        if (buf.readInt8() == 1) {
            symbols = ByteData.toFloat32List(Zlib.inflate(buf.readUint8List()));
            repeats = readRepeatsArchived(buf);
        } else {
            symbols = buf.readFloat32List();
            repeats = readRepeats(buf);
        }
        int idxDecoded = 0;
        int idxSymbols = 0;
        for (int n = 0; n < repeats.length; n++) {
            int cnt = repeats[n] & (~BLOCK_FLAG_MASK);
            boolean isBlock = (repeats[n] & BLOCK_FLAG_MASK) == BLOCK_FLAG_MASK;
            if (isBlock) {
                while (cnt-- > 0)
                    decoded[idxDecoded++] = symbols[idxSymbols++];
            } else {
                while (cnt-- > 0)
                    decoded[idxDecoded++] = symbols[idxSymbols];
                idxSymbols++;
            }
        }
        return decoded;
    }

    public static double[] decode64(BufferAccessor buf) {
        double[] decoded = new double[buf.readInt32()];
        double[] symbols;
        int[] repeats;
        if (buf.readInt8() == 1) {
            symbols = ByteData.toFloat64List(Zlib.inflate(buf.readUint8List()));
            repeats = readRepeatsArchived(buf);
        } else {
            symbols = buf.readFloat64List();
            repeats = readRepeats(buf);
        }
        int idxDecoded = 0;
        int idxSymbols = 0;
        for (int n = 0; n < repeats.length; n++) {
            int cnt = repeats[n] & (~BLOCK_FLAG_MASK);
            boolean isBlock = (repeats[n] & BLOCK_FLAG_MASK) == BLOCK_FLAG_MASK;
            if (isBlock) {
                while (cnt-- > 0)
                    decoded[idxDecoded++] = symbols[idxSymbols++];
            } else {
                while (cnt-- > 0)
                    decoded[idxDecoded++] = symbols[idxSymbols];
                idxSymbols++;
            }
        }
        return decoded;
    }

    // Repeats are a uint8 list; return them as unsigned ints (0..255).
    private static int[] readRepeats(BufferAccessor buf) {
        return toUnsigned(buf.readUint8List());
    }

    private static int[] readRepeatsArchived(BufferAccessor buf) {
        return toUnsigned(Zlib.inflate(buf.readUint8List()));
    }

    private static int[] toUnsigned(byte[] b) {
        int[] r = new int[b.length];
        for (int i = 0; i < b.length; i++)
            r[i] = b[i] & 0xFF;
        return r;
    }
}
