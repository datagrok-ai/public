package serialization;

import java.util.*;
import java.util.regex.*;
import java.util.function.*;
import java.util.stream.Collectors;


// Data frame.
public class DataFrame {
    public String name;
    public Integer rowCount = 0;
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
        String collect = columns.stream()
                .map(column -> column.name)
                .collect(Collectors.joining(delimiter));
        buffer.append(collect);
        buffer.append(System.lineSeparator());

        for (int r = 0; r < rowCount; r++) {
            final int index = r;
            String row = columns.stream()
                    .map(column -> columnToStr(column, index))
                    .collect(Collectors.joining(delimiter));
            buffer.append(row);
            buffer.append(System.lineSeparator());
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

    public void merge(DataFrame dataFrame) {
        if (!Objects.equals(columns.size(), dataFrame.columns.size()) && rowCount != 0) {
            throw new RuntimeException("Can't merge dataframes with different row count");
        }
        if (rowCount == 0) {
            addColumns(dataFrame.columns);
            return;
        }
        for (int i = 0; i < columns.size(); i++) {
            Column column = columns.get(i);
            Column columnToMerge = dataFrame.columns.get(i);
            if (!column.name.equals(columnToMerge.name)) {
                throw new RuntimeException("Can't merge dataframes of columns with different names");
            }
            for (int j = 0; j < columnToMerge.length; j++) {
                column.add(columnToMerge.get(j));
            }
        }
        Column column = columns.get(0);
        if (column != null) {
            rowCount = column.length;
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        DataFrame dataFrame = (DataFrame) o;
        return Objects.equals(name, dataFrame.name) && Objects.equals(rowCount, dataFrame.rowCount) && Objects.equals(columns, dataFrame.columns) && Objects.equals(tags, dataFrame.tags);
    }

    @Override
    public int hashCode() {
        return Objects.hash(name, rowCount, columns, tags);
    }

    public int memoryInBytes() {
        int sum = 0;
        for (Column col : columns)
            sum += col.memoryInBytes();
        return sum;
    }
}
