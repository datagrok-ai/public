package serialization.codecs;

import java.nio.charset.StandardCharsets;

import serialization.BufferAccessor;
import serialization.Column;

// Decode-only port of StringsSquash.decode
// (ddt/lib/src/serialization/string_squash.dart:70-90).
//
// The strings all share the same length. The first string is the reference; for
// each character position that differs across rows (compare[n] != 0), a nested
// int column holds each row's code unit at that position, UTF-8-decoded to a
// string whose n-th char replaces position n of row n. Note the Dart quirk that
// the per-position code units are UTF-8-decoded as bytes (faithful for ASCII).
public class StringSquash {
    public static String[] decode(BufferAccessor buf) {
        int count = buf.readInt32();
        String first = buf.readString();
        byte[] compare = buf.readUint8List();

        char[][] rows = new char[count][];
        for (int m = 0; m < count; m++)
            rows[m] = first.toCharArray();

        for (int n = 0; n < compare.length; n++) {
            if ((compare[n] & 0xFF) != 0) {
                Column<?> col = buf.readColumn();
                int len = col.getLength();
                int[] codeUnits = (int[]) col.toArray();
                byte[] bytes = new byte[len];
                for (int m = 0; m < len; m++)
                    bytes[m] = (byte) codeUnits[m];
                String str = new String(bytes, StandardCharsets.UTF_8);
                for (int m = 1; m < count; m++)
                    rows[m][n] = str.charAt(m);
            }
        }

        String[] out = new String[count];
        for (int m = 0; m < count; m++)
            out[m] = new String(rows[m]);
        return out;
    }
}
