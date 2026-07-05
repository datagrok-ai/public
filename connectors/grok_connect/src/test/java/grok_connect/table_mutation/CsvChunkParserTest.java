package grok_connect.table_mutation;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * Unit tests for the streaming RFC 4180 tokenizer (connector-writes WO-5): chunk-boundary splitting
 * of quoted fields carrying commas/newlines/quotes, CRLF handling, the null-vs-empty-string
 * distinction, and multi-byte UTF-8 characters split across a chunk boundary.
 */
class CsvChunkParserTest {
    /** Feeds each chunk in turn (order matters — this is the whole point), then finishes. */
    private List<List<String>> parse(byte[]... chunks) {
        CsvChunkParser parser = new CsvChunkParser();
        List<List<String>> rows = new ArrayList<>();
        for (byte[] chunk : chunks)
            rows.addAll(parser.feed(chunk));
        rows.addAll(parser.finish());
        return rows;
    }

    private byte[] b(String s) {
        return s.getBytes(StandardCharsets.UTF_8);
    }

    @DisplayName("Simple record and finish flushes the last row without a trailing newline")
    @Test
    public void simpleRecords() {
        List<List<String>> rows = parse(b("a,b,c\nd,e,f"));
        Assertions.assertEquals(2, rows.size());
        Assertions.assertEquals(java.util.Arrays.asList("a", "b", "c"), rows.get(0));
        Assertions.assertEquals(java.util.Arrays.asList("d", "e", "f"), rows.get(1));
    }

    @DisplayName("Empty unquoted field is null; empty quoted field is the empty string")
    @Test
    public void nullVsEmptyString() {
        List<List<String>> rows = parse(b("a,,\"\",d\n"));
        List<String> row = rows.get(0);
        Assertions.assertEquals("a", row.get(0));
        Assertions.assertNull(row.get(1));          // empty unquoted -> null
        Assertions.assertEquals("", row.get(2));    // "" -> empty string
        Assertions.assertEquals("d", row.get(3));
    }

    @DisplayName("A quoted field containing a comma and newline, split across a chunk boundary, reassembles")
    @Test
    public void quotedFieldWithCommaNewlineSplitAcrossChunks() {
        // record: "x,\ny",z  — boundary falls inside the quoted field, right after the comma
        List<List<String>> rows = parse(b("\"x,"), b("\ny\",z\n"));
        Assertions.assertEquals(1, rows.size());
        Assertions.assertEquals("x,\ny", rows.get(0).get(0));
        Assertions.assertEquals("z", rows.get(0).get(1));
    }

    @DisplayName("Escaped double-quotes, with the boundary landing between the escape pair")
    @Test
    public void escapedQuotesSplitAcrossChunks() {
        // record: "a""b",z  split as `"a"` + `"b",z\n` (boundary between the two quotes of the escape)
        List<List<String>> rows = parse(b("\"a\""), b("\"b\",z\n"));
        Assertions.assertEquals("a\"b", rows.get(0).get(0));
        Assertions.assertEquals("z", rows.get(0).get(1));
    }

    @DisplayName("CRLF terminators, including a boundary between CR and LF")
    @Test
    public void crlf() {
        Assertions.assertEquals(2, parse(b("a,b\r\nc,d\r\n")).size());
        List<List<String>> split = parse(b("a,b\r"), b("\nc,d\n"));
        Assertions.assertEquals(2, split.size());
        Assertions.assertEquals(java.util.Arrays.asList("a", "b"), split.get(0));
        Assertions.assertEquals(java.util.Arrays.asList("c", "d"), split.get(1));
    }

    @DisplayName("A multi-byte UTF-8 character split across a chunk boundary is decoded intact")
    @Test
    public void multiByteCharSplitAcrossChunks() {
        byte[] full = b("café,→\n"); // é = 2 bytes, → = 3 bytes
        // split so that the é (bytes 3-4) straddles the boundary
        byte[] first = new byte[4];
        byte[] second = new byte[full.length - 4];
        System.arraycopy(full, 0, first, 0, 4);
        System.arraycopy(full, 4, second, 0, second.length);
        List<List<String>> rows = parse(first, second);
        Assertions.assertEquals("café", rows.get(0).get(0));
        Assertions.assertEquals("→", rows.get(0).get(1));
    }

    @DisplayName("Blank lines between records are skipped")
    @Test
    public void blankLinesSkipped() {
        List<List<String>> rows = parse(b("a,b\n\nc,d\n"));
        Assertions.assertEquals(2, rows.size());
    }

    @DisplayName("A newline embedded in a quoted field does not terminate the record")
    @Test
    public void embeddedNewlineDoesNotSplitRecord() {
        List<List<String>> rows = parse(b("\"line1\nline2\",x\n"));
        Assertions.assertEquals(1, rows.size());
        Assertions.assertEquals("line1\nline2", rows.get(0).get(0));
    }

    @DisplayName("Byte-by-byte feeding produces the same result as a single feed")
    @Test
    public void byteByByteFeeding() {
        String csv = "\"a,b\",\"c\"\"d\",\ne,f\n";
        CsvChunkParser oneShot = new CsvChunkParser();
        List<List<String>> expected = new ArrayList<>(oneShot.feed(b(csv)));
        expected.addAll(oneShot.finish());

        CsvChunkParser drip = new CsvChunkParser();
        List<List<String>> actual = new ArrayList<>();
        for (byte by : b(csv)) {
            ByteArrayOutputStream one = new ByteArrayOutputStream();
            one.write(by);
            actual.addAll(drip.feed(one.toByteArray()));
        }
        actual.addAll(drip.finish());
        Assertions.assertEquals(expected, actual);
    }
}
