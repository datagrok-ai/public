package serialization;

import java.io.ByteArrayOutputStream;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

// Zlib (deflate) inflate helper for the archive-flag decode branches.
// Mirrors Dart's zip.ZLibDecoder().decodeBytes(). Current writers force the
// archive flag to 0, so this is only exercised by legacy blobs.
public class Zlib {
    public static byte[] inflate(byte[] data) {
        Inflater inflater = new Inflater();
        try {
            inflater.setInput(data);
            ByteArrayOutputStream out = new ByteArrayOutputStream(Math.max(64, data.length * 2));
            byte[] chunk = new byte[4096];
            while (!inflater.finished()) {
                int n = inflater.inflate(chunk);
                if (n == 0) {
                    if (inflater.finished() || inflater.needsDictionary() || inflater.needsInput())
                        break;
                }
                out.write(chunk, 0, n);
            }
            return out.toByteArray();
        } catch (DataFormatException e) {
            throw new RuntimeException("Zlib inflate failed", e);
        } finally {
            inflater.end();
        }
    }
}
