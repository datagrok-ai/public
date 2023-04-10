package serialization;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DataFrame {
    private static final String delimiter = ",";
    public String name;
    public Integer rowCount = 0;
    public List<Column> columns = new ArrayList<>();
    public Map<String, String> tags;

    public DataFrame() {
    }

    public void addColumn(Column col) {
        rowCount = col.length;
        col.name = getUniqueName(col.name);
        columns.add(col);
    }

    public void addColumns(List<Column> cols) {
        if (cols.size() > 0) {
            Column column = cols.get(0);
            if (column.getType().equals(Types.COLUMN_LIST)) {
                Column column1 = (Column) column.get(0);
                rowCount = column1 == null ? 0 : (column1).length;
            } else {
                rowCount = column.length;
            }
            addColumnsRecursive(cols);
        }
    }

    private void addColumnsRecursive(List<Column> cols) {
        for (Column col : cols) {
            if (col == null) break;
            if (col.getType().equals(Types.COLUMN_LIST)) {
                ComplexTypeColumn complexTypeColumn = (ComplexTypeColumn) col;
                addColumnsRecursive(new ArrayList<>(Arrays.asList(complexTypeColumn.getAll())));
            } else {
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
