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

    public byte[] toByteArray() {
        DataFrame[] tables = {this};
        return (new TablesBlob(tables)).toByteArray();
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
