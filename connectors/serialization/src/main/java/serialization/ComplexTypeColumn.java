package serialization;

import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class ComplexTypeColumn extends Column<List<Column>> {
    private static final String DATA_FIELD_NAME = "data";
    private static final int DEFAULT_ARRAY_SIZE = 15;
    private Column[] data;

    public ComplexTypeColumn() {
        this.data = new Column[DEFAULT_ARRAY_SIZE];
    }

    public ComplexTypeColumn(Column[] columns) {
        this.data = columns;
    }

    @Override
    public String getType() {
        return Types.COLUMN_LIST;
    }

    @Override
    public void encode(BufferAccessor buf) {
        for (int i = 0; i < length; i++) {
            Column column = data[i];
            if (column != null) {
                column.encode(buf);
            }
        }
    }

    @Override
    public void add(List<Column> columns) {
        ensureSpace(columns.size());
        if (length == 0) {
            System.arraycopy(columns.toArray(new Column[0]), 0, data, 0, columns.size());
            length += columns.size();
        } else {
            for (Column column: columns) {
                String name = column.name;
                Optional<Column> first = Arrays.stream(data)
                        .filter(column1 -> column1.name.equals(name))
                        .findFirst();
                if (first.isPresent()) {
                    Column column1 = first.get();
                    if (column1.getType().equals(Types.BOOL)) {
                        int value = (int) column.get(0);
                        column1.add(value == 1);
                    } else {
                        column1.add(column.get(0));
                    }
                } else {
                    int lengthDiff = Math.abs(data[0].length - column.length);
                    if (lengthDiff > 0) {
                        insertNull(column, lengthDiff);
                        appendNull(data, lengthDiff);
                    }
                    data[length++] = column;
                }
            }
        }
    }

    @Override
    public void addAll(List<Column>[] value) {
        throw new UnsupportedOperationException("Not supported");
    }


    @Override
    public Object get(int idx) {
        if (idx >= data.length) {
            return null;
        }
        return data[idx];
    }

    @Override
    public long memoryInBytes() {
        if (data[0] == null) {
            return 0L;
        }
        return Arrays.stream(data)
                .map(Column::memoryInBytes)
                .reduce(0L, Long::sum);
    }

    @Override
    public boolean isNone(int idx) {
        return idx >= data.length;
    }

    @Override
    public void empty() {
        data = new Column[DEFAULT_ARRAY_SIZE];
    }

    public Column[] getAll() {
        Column[] returnArray = new Column[length];
        System.arraycopy(data, 0, returnArray, 0, length);
        return returnArray;
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            Column[] newData = new Column[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    private void insertNull(Column column, int nullCount) {
        try {
            Field data = column.getClass().getDeclaredField(DATA_FIELD_NAME);
            data.setAccessible(true);
            Class<?> componentType = data.getClass().getComponentType();
            Object array = Array.newInstance(componentType, column.length + nullCount);
            for (int i = 0; i < nullCount; i++) {
                Array.set(array, i, null);
            }
            for (int j = nullCount, i = 0; i < column.length; j++, i++) {
                Array.set(array, j, column.get(i));
            }
            data.set(column, array);
            data.setAccessible(false);
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new RuntimeException("Something went wrong when inserting null", e);
        }
    }

    private void insertNull(Column[] columns, int nullCount) {
        for (Column column: columns) {
            insertNull(column, nullCount);
        }
    }

    private void appendNull(Column column, int nullCount) {
        for (int i = 0; i < nullCount; i++) {
            column.add(null);
        }
    }

    private void appendNull(Column[] columns, int nullCount) {
        for (Column column: columns) {
            appendNull(column, nullCount);
        }
    }
}
