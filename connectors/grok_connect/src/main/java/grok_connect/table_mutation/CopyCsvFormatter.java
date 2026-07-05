package grok_connect.table_mutation;

import java.nio.charset.StandardCharsets;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.List;

import serialization.Column;
import serialization.DataFrame;
import serialization.DateTimeColumn;
import serialization.FloatColumn;
import serialization.Types;

/**
 * Formats a decoded d42 chunk into {@code COPY ... FROM STDIN (FORMAT csv)} text (java-d42-reader WO-4):
 * a none cell → an empty unquoted field (COPY's NULL), strings RFC-4180-quoted only when they carry a
 * comma/quote/newline, {@code bool} → {@code true|false}, {@code double} via {@link Double#toString}
 * (COPY accepts scientific notation), {@code datetime} → {@code yyyy-MM-dd HH:mm:ss.SSSSSS} UTC from
 * epoch-µs. Column order follows the chunk (which the manager validated against the header order), so
 * the emitted fields align with the {@code COPY} column list. The typed binder ({@link MutationRunner})
 * and this formatter must store identical values for the same chunk (loader-parity pin).
 */
class CopyCsvFormatter {
    private static final DateTimeFormatter TIMESTAMP =
            DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss.SSSSSS");

    /** Formats every row of {@code chunk} as UTF-8 COPY-csv text (LF-terminated rows, no header). */
    static byte[] format(DataFrame chunk) {
        List<Column<?>> columns = chunk.getColumns();
        StringBuilder sb = new StringBuilder(chunk.rowCount * Math.max(1, columns.size()) * 8);
        for (int r = 0; r < chunk.rowCount; r++) {
            for (int c = 0; c < columns.size(); c++) {
                if (c > 0)
                    sb.append(',');
                appendCell(sb, columns.get(c), r);
            }
            sb.append('\n');
        }
        return sb.toString().getBytes(StandardCharsets.UTF_8);
    }

    private static void appendCell(StringBuilder sb, Column<?> col, int row) {
        if (col.isNone(row))
            return; // empty unquoted field = COPY csv NULL
        switch (col.getType()) {
            case Types.INT:
            case Types.BIG_INT:
                sb.append(col.get(row).toString());
                break;
            case Types.NUM:
            case Types.FLOAT:
                sb.append(Double.toString(((FloatColumn) col).getDouble(row)));
                break;
            case Types.BOOL:
                sb.append(((Boolean) col.get(row)) ? "true" : "false");
                break;
            case Types.DATE_TIME:
                sb.append(formatDateTime(((DateTimeColumn) col).get(row)));
                break;
            default:
                sb.append(quote(col.get(row).toString()));
        }
    }

    private static String formatDateTime(double micros) {
        long us = (long) micros;
        long seconds = Math.floorDiv(us, 1_000_000L);
        long microOfSecond = Math.floorMod(us, 1_000_000L);
        return LocalDateTime.ofEpochSecond(seconds, (int) (microOfSecond * 1000L), ZoneOffset.UTC).format(TIMESTAMP);
    }

    private static String quote(String s) {
        if (s.indexOf(',') < 0 && s.indexOf('"') < 0 && s.indexOf('\n') < 0 && s.indexOf('\r') < 0)
            return s;
        return "\"" + s.replace("\"", "\"\"") + "\"";
    }
}
