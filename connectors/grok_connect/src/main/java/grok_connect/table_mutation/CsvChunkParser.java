package grok_connect.table_mutation;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * Streaming RFC 4180 CSV tokenizer (UTF-8, no header row) for the bulk mutation transport.
 * Fed one byte chunk at a time; emits only the records completed so far and carries the partial
 * last record across chunk boundaries — so a quoted field containing commas or newlines split by a
 * chunk boundary reassembles correctly. Structural characters (comma, quote, CR, LF) are single
 * bytes in UTF-8 and never collide with multi-byte sequences, so scanning at the byte level keeps
 * multi-byte characters intact even when a chunk boundary splits them.
 *
 * <p>Null convention (matches Postgres {@code COPY ... FORMAT csv}): an empty <em>unquoted</em> field
 * is {@code null}; an empty <em>quoted</em> field ({@code ""}) is the empty string.
 */
class CsvChunkParser {
    private final List<String> fields = new ArrayList<>();
    private final ByteArrayOutputStream field = new ByteArrayOutputStream();
    private boolean quoted;       // current field opened with a quote
    private boolean hasBytes;     // current field received any bytes
    private boolean dirty;        // current record has any content (distinguishes a blank line)
    private boolean inQuotes;
    private boolean quotePending; // inside quotes, saw a '"' — next byte decides escape vs close
    private boolean crPending;    // saw CR outside quotes — consume a following LF as part of CRLF

    /** Parses whatever records the accumulated bytes now complete; retains the partial tail. */
    List<List<String>> feed(byte[] chunk) {
        List<List<String>> out = new ArrayList<>();
        for (byte b : chunk)
            process(b, out);
        return out;
    }

    /** Flushes the final record if the stream ended without a trailing newline. */
    List<List<String>> finish() {
        List<List<String>> out = new ArrayList<>();
        if (crPending) {
            crPending = false;
            endRecord(out);
        }
        else if (quotePending) {
            quotePending = false;
            inQuotes = false;
            endRecord(out);
        }
        else if (dirty || hasBytes)
            endRecord(out);
        return out;
    }

    private void process(byte b, List<List<String>> out) {
        if (crPending) {
            crPending = false;
            endRecord(out);
            if (b == '\n')
                return; // swallow the LF of a CRLF terminator
        }
        if (quotePending) {
            quotePending = false;
            if (b == '"') {
                field.write('"'); // escaped quote inside a quoted field
                return;
            }
            inQuotes = false; // the pending quote closed the field; reprocess b unquoted
        }
        if (inQuotes) {
            if (b == '"')
                quotePending = true;
            else
                field.write(b);
            return;
        }
        switch (b) {
            case '"':
                quoted = true;
                inQuotes = true;
                dirty = true;
                break;
            case ',':
                endField();
                dirty = true;
                break;
            case '\r':
                crPending = true;
                break;
            case '\n':
                endRecord(out);
                break;
            default:
                field.write(b);
                hasBytes = true;
                dirty = true;
        }
    }

    private void endField() {
        String value = quoted || hasBytes ? new String(field.toByteArray(), StandardCharsets.UTF_8) : null;
        fields.add(value);
        field.reset();
        quoted = false;
        hasBytes = false;
    }

    private void endRecord(List<List<String>> out) {
        if (!dirty) { // blank line — no content
            reset();
            return;
        }
        endField();
        out.add(new ArrayList<>(fields));
        reset();
    }

    private void reset() {
        fields.clear();
        field.reset();
        quoted = false;
        hasBytes = false;
        dirty = false;
    }
}
