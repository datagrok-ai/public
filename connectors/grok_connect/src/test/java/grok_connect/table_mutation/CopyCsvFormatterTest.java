package grok_connect.table_mutation;

import java.nio.charset.StandardCharsets;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import serialization.BigIntColumn;
import serialization.BoolColumn;
import serialization.Column;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;

/**
 * Unit tests for {@link CopyCsvFormatter} (java-d42-reader WO-4): typed d42 columns → {@code COPY csv}
 * text. No Docker — asserts the exact bytes the Postgres COPY fast path streams.
 */
public class CopyCsvFormatterTest {
    private String[] lines(Column<?>... cols) {
        String text = new String(CopyCsvFormatter.format(DataFrame.fromColumns(cols)), StandardCharsets.UTF_8);
        // format ends every row with '\n'; drop the trailing empty element from split
        return text.split("\n", -1);
    }

    /** Formats a single-column single-cell frame and returns the field text (a quoted field may hold a '\n'). */
    private String oneCell(Column<?> col) {
        String text = new String(CopyCsvFormatter.format(DataFrame.fromColumns(col)), StandardCharsets.UTF_8);
        return text.substring(0, text.length() - 1); // strip the trailing row terminator
    }

    @DisplayName("none cells become empty unquoted fields (COPY csv NULL)")
    @Test
    public void noneIsEmptyUnquoted() {
        String[] l = lines(new IntColumn("i", new Integer[] {1, null, 3}));
        Assertions.assertEquals("1", l[0]);
        Assertions.assertEquals("", l[1]); // IntColumn.None -> empty
        Assertions.assertEquals("3", l[2]);
        // a null (empty '') string category is none too
        String[] s = lines(new StringColumn("s", new String[] {"x", null, ""}));
        Assertions.assertEquals("x", s[0]);
        Assertions.assertEquals("", s[1]);
        Assertions.assertEquals("", s[2]); // '' is the string None slot
    }

    @DisplayName("RFC-4180 quoting only when the string carries comma / quote / newline")
    @Test
    public void rfc4180Quoting() {
        Assertions.assertEquals("plain", oneCell(new StringColumn("s", new String[] {"plain"})));
        Assertions.assertEquals("\"a,b\"", oneCell(new StringColumn("s", new String[] {"a,b"})));
        Assertions.assertEquals("\"he said \"\"hi\"\"\"", oneCell(new StringColumn("s", new String[] {"he said \"hi\""})));
        Assertions.assertEquals("\"line1\nline2\"", oneCell(new StringColumn("s", new String[] {"line1\nline2"})));
        Assertions.assertEquals("\"a\rb\"", oneCell(new StringColumn("s", new String[] {"a\rb"})));
        Assertions.assertEquals("tab\tok", oneCell(new StringColumn("s", new String[] {"tab\tok"}))); // tab is not structural
    }

    @DisplayName("bool -> true|false; int/bigint verbatim; double via Double.toString incl. NaN/-0.0/scientific")
    @Test
    public void scalarFormatting() {
        String[] b = lines(new BoolColumn("b", new Boolean[] {true, false}));
        Assertions.assertEquals("true", b[0]);
        Assertions.assertEquals("false", b[1]);

        String[] big = lines(new BigIntColumn("g", new String[] {"9223372036854775807", "-42"}));
        Assertions.assertEquals("9223372036854775807", big[0]);
        Assertions.assertEquals("-42", big[1]);

        String[] d = lines(FloatColumn.double64("d", new Double[] {
                1.5, 1.0000000000000002, Double.NaN, -0.0, 1e308}));
        Assertions.assertEquals("1.5", d[0]);
        Assertions.assertEquals("1.0000000000000002", d[1]);
        Assertions.assertEquals("NaN", d[2]);
        Assertions.assertEquals("-0.0", d[3]);
        Assertions.assertEquals(Double.toString(1e308), d[4]);
    }

    @DisplayName("datetime -> 'yyyy-MM-dd HH:mm:ss.SSSSSS' UTC from epoch-us (us precision, pre-1970)")
    @Test
    public void dateTimeFormatting() {
        String[] l = lines(new DateTimeColumn("t", new Double[] {
                1_000_000d,          // 1970-01-01 00:00:01.000000
                1_500_000d,          // 1970-01-01 00:00:01.500000
                1_000_123_456d,      // 1970-01-01 00:16:40.123456 (us precision)
                -1_000_000d,         // 1969-12-31 23:59:59.000000
                -500_000d}));        // 1969-12-31 23:59:59.500000
        Assertions.assertEquals("1970-01-01 00:00:01.000000", l[0]);
        Assertions.assertEquals("1970-01-01 00:00:01.500000", l[1]);
        Assertions.assertEquals("1970-01-01 00:16:40.123456", l[2]);
        Assertions.assertEquals("1969-12-31 23:59:59.000000", l[3]);
        Assertions.assertEquals("1969-12-31 23:59:59.500000", l[4]);
    }

    @DisplayName("multi-column rows assemble comma-joined in column order, nulls interspersed")
    @Test
    public void multiColumnRowAssembly() {
        String[] l = lines(
                new IntColumn("id", new Integer[] {1, 2}),
                new StringColumn("note", new String[] {"a,b", null}),
                new BoolColumn("active", new Boolean[] {true, false}));
        Assertions.assertEquals("1,\"a,b\",true", l[0]);
        Assertions.assertEquals("2,,false", l[1]);
    }
}
