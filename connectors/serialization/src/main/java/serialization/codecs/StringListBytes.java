package serialization.codecs;

import java.nio.charset.StandardCharsets;

// Decode-only port of _StringListToUint8List.decode
// (ddt/lib/src/data_frame/encoders/string_column_encoders.dart:294-313).
// Used by the string:prefixes and string:zlib decoders.
//
// Two layouts:
//   * transpose - all strings share one UTF-8 byte length (lengths[0]); bytes are
//     stored column-major: byte m of string n is at bytes[m * length + n].
//   * flat - each string's UTF-8 byte length is lengths[n]; bytes are concatenated.
public class StringListBytes {
    public static String[] decode(byte[] bytes, int length, int[] lengths, boolean transpose) {
        String[] strings = new String[length];
        if (transpose) {
            byte[] stringBytes = new byte[lengths[0]];
            for (int n = 0; n < length; n++) {
                for (int m = 0; m < lengths[0]; m++)
                    stringBytes[m] = bytes[m * length + n];
                strings[n] = new String(stringBytes, StandardCharsets.UTF_8);
            }
        } else {
            int offset = 0;
            for (int n = 0; n < length; n++) {
                strings[n] = new String(bytes, offset, lengths[n], StandardCharsets.UTF_8);
                offset += lengths[n];
            }
        }
        return strings;
    }
}
