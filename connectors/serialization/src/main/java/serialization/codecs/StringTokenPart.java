package serialization.codecs;

import serialization.IntColumn;
import serialization.StringColumn;

// Port of StringTokenPart (ddt/lib/src/serialization/string_tokens.dart, GROK-20761);
// same structure and names. One part of a StringTokens template: a literal separator,
// an int slot, a long-int slot, or a string slot.
public class StringTokenPart {
    public static final int LITERAL = 0;
    public static final int INT = 1;
    public static final int STRING = 2;
    public static final int LONG_INT = 3;

    public static final int FLAG_PADDED = 1;
    public static final int FLAG_HEX_UPPER = 2;
    public static final int FLAG_HEX_LOWER = 4;

    public int kind;
    public String text;
    public int flags = 0;
    public int width = 0;

    // Extraction coordinates: class run within the template field (-1 = whole field) and,
    // for digit runs split into several parts, the offset of this part inside the run.
    public int run = -1;
    public int chunkOffset = -1;

    public int[] ints;
    public long[] longs;
    public String[] strings;
    public IntColumn intColumn;       // decoded INT payload / estimate-phase writer column
    public StringColumn stringColumn; // decoded STRING payload / estimate-phase writer column

    // long parts only: 0 = first + bit-packed deltas, 1 = hi/lo int columns.
    public int longMode;
    public long[] longDeltaRange;
    public int[] hi, lo;

    public StringTokenPart(int kind) {
        this.kind = kind;
    }

    public boolean isHex() {
        return (flags & (FLAG_HEX_UPPER | FLAG_HEX_LOWER)) != 0;
    }

    @Override public String toString() {
        return kind == LITERAL ? text
            : kind == STRING ? "{str}"
            : kind == LONG_INT ? "{long:" + width + "}"
            : isHex() ? "{hex:" + width + "}"
            : (flags & FLAG_PADDED) != 0 ? "{int:" + width + "}" : "{int}";
    }
}
