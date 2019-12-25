package serialization;

import java.util.*;


// Data frame.
public class DataFrame {
    public String name;
    public Integer rowCount;
    public List<Column> columns = new ArrayList<>();
    public Map<String, String> tags;

    public DataFrame() {
    }

    public void addColumn(Column col) {
        rowCount = col.length;
        columns.add(col);
    }

    public void addColumns(List<Column> cols) {
        if (cols.size() > 0) {
            rowCount = cols.get(0).length;
            columns.addAll(cols);
        }
    }

    public byte[] toByteArray() {
        DataFrame[] tables = {this};
        return (new TablesBlob(tables)).toByteArray();
    }
}
