package serialization;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class DataFrame implements Taggable {
    private static final String delimiter = ",";
    private final List<Column<?>> columns = new ArrayList<>();
    private final Map<String, String> tags = new HashMap<>();

    public String name;
    public Integer rowCount = 0;

    public DataFrame() {
    }

    public static DataFrame fromColumns(Column<?>...columns) {
        DataFrame df = new DataFrame();
        df.addColumns(columns);
        return df;
    }

    public List<Column<?>> getColumns() {
        return columns;
    }

    public Column<?> getColumn(String name) {
        for (Column<?> col : columns)
            if (col.getName().equals(name))
                return col;
        return null;
    }

    public Column<?> getColumn(int index) {
        return columns.get(index);
    }

    public int getColumnCount() {
        return columns.size();
    }

    public void removeColumn(String name) {
        columns.removeIf(c -> c.getName().equalsIgnoreCase(name));
    }

    @SuppressWarnings("unchecked")
    public void addRow(Object... objects) {
        if (objects.length != columns.size())
            throw new RuntimeException("Objects length does not match columns length");
        for (int i = 0; i < columns.size(); i++)
            ((Column<Object>) columns.get(i)).add(objects[i]);
        rowCount++;
    }

    public void addColumn(Column<?> col) {
        rowCount = col.getLength();
        col.setName(getUniqueName(col.getName()));
        columns.add(col);
    }

    public void addColumns(Column<?>[] cols) {
        if (cols != null && cols.length > 0) {
            Column<?> column = cols[0];
            if (column.getType().equals(Types.COLUMN_LIST)) {
                ComplexTypeColumn ctc = (ComplexTypeColumn) column;
                Column<?>[] inner = ctc.convertToColumns();
                rowCount = inner.length > 0 ? inner[0].getLength() : 0;
            } else {
                rowCount = column.getLength();
            }
            addColumnsRecursive(cols);
        }
    }

    private void addColumnsRecursive(Column<?>[] cols) {
        for (Column<?> col : cols) {
            if (col == null) break;
            if (col.getType().equals(Types.COLUMN_LIST)) {
                ComplexTypeColumn complexTypeColumn = (ComplexTypeColumn) col;
                addColumnsRecursive(complexTypeColumn.convertToColumns());
            } else {
                col.setName(getUniqueName(col.getName()));
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
                .map(Column::getName)
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

    private String columnToStr(Column<?> col, int row) {
        if (col.isNone(row))
            return "";
        String str = String.valueOf(col.get(row));
        return csvQuote(str);
    }

    private String getUniqueName(String name) {
        Predicate<String> unique = (candidate) -> {
            for (Column<?> column : columns)
                if (column.getName().equalsIgnoreCase(candidate))
                    return false;
            return true;
        };

        if (unique.test(name))
            return name;

        int i = 1;
        String core = name;

        Pattern pattern = Pattern.compile("^(.*) \\((\\d+)\\)$");
        Matcher matcher = pattern.matcher(name);
        if (matcher.find()) {
            core = matcher.group(1);
            i = Integer.parseInt(matcher.group(2));
        }

        for (; true; i++) {
            String newName = core + " (" + i + ")";
            if (unique.test(newName))
                return newName;
        }
    }

    @SuppressWarnings("unchecked")
    public void merge(DataFrame dataFrame) {
        if (!Objects.equals(columns.size(), dataFrame.columns.size()) && rowCount != 0)
            throw new RuntimeException("Can't merge dataframes with different row count");
        if (rowCount == 0) {
            addColumns(dataFrame.columns.toArray(new Column<?>[0]));
            return;
        }
        for (int i = 0; i < columns.size(); i++) {
            Column<Object> column = (Column<Object>) columns.get(i);
            Column<?> columnToMerge = dataFrame.columns.get(i);
            if (!column.getName().equals(columnToMerge.getName()))
                throw new RuntimeException("Can't merge dataframes of columns with different names");
            for (int j = 0; j < columnToMerge.getLength(); j++)
                column.add(columnToMerge.get(j));
        }
        Column<?> column = columns.get(0);
        if (column != null)
            rowCount = column.getLength();
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

    @Override
    public Map<String, String> getTags() {
        return tags;
    }

    public long memoryInBytes() {
        long sum = 0;
        for (Column<?> col : columns)
            sum += col.memoryInBytes();
        return sum;
    }
}
