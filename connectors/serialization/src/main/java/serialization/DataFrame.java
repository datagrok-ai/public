package serialization;

import java.util.*;
import java.util.regex.*;
import java.util.function.*;


// Data frame.
public class DataFrame {
    public String name;
    public Integer rowCount;
    public List<Column> columns = new ArrayList<>();
    public Map<String, String> tags;

    private static final String delimiter = ",";

    public DataFrame() {
    }

    public void addColumn(Column col) {
        rowCount = col.length;
        col.name = getUniqueName(col.name);
        columns.add(col);
    }

    public void addColumns(List<Column> cols) {
        if (cols.size() > 0) {
            rowCount = cols.get(0).length;
            for (Column col : cols) {
                col.name = getUniqueName(col.name);
                columns.add(col);
            }
        }
    }

    String csvQuote(String s) {
        boolean quotes = s.contains("\"");
        s = s.replaceAll("\"", "\"\"");
        return quotes || s.contains(delimiter) || s.contains("\n") ? "\"" + s + "\"" : s;
    }

    public String toCsv() {
        StringBuilder buffer = new StringBuilder();

        for (Column column : columns) {
            buffer.append(column.name);
            buffer.append(delimiter);
        }
        buffer.append("\n");

        for (int r = 0; r < columns.get(0).length; r++) {
            for (Column column : columns) {
                buffer.append(columnToStr(column, r));
                buffer.append(delimiter);
            }
            buffer.append("\n");
        }

        return buffer.toString();
    }

    public byte[] toByteArray() {
        DataFrame[] tables = {this};
        return (new TablesBlob(tables)).toByteArray();
    }

    private String columnToStr(Column col, int row) {
        if (col.isNone(row))
            return "";
        String str = String.valueOf(col.get(row));
        return csvQuote(str);
    }

    private String getUniqueName(String name) {
        Predicate<String> unique = (candidate) -> {
            String lowerCandidate = candidate.toLowerCase();
            for (Column column : columns)
                if (column.name.equals(lowerCandidate))
                    return false;
            return true;
        };

        if (unique.test(name))
            return name;

        int i = 1;
        String core = name;

        Pattern pattern = Pattern.compile("^(.*) \\(([0-9]+)\\)$");
        Matcher matcher = pattern.matcher(name);
        if (matcher.find()) {
            core = matcher.group(1);
            i = Integer.parseInt(matcher.group(2));
        }

        for (; true; i++) {
            String newName = core + " (" + String.valueOf(i) + ")";
            if (unique.test(newName))
                return newName;
        }
    }
}
