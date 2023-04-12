package serialization;

import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

public class ComplexTypeColumn extends Column<List<Column>> {
    private static final String DATA_FIELD_NAME = "data";
    private static final int DEFAULT_ARRAY_SIZE = 100;
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
        processListMap(columns);
        if (length == 0) {
            System.arraycopy(columns.toArray(new Column[0]), 0, data, 0, columns.size());
            length += columns.size();
        } else {
            // append nulls for each column in case there are no such column
            appendNull(data, 1);
            for (Column column: columns) {
                String name = column.name;
                Column first = null;
                for (Column col: data) {
                    if (col != null && col.name.equals(name)) {
                        first = col;
                        break;
                    }
                }
                if (first != null) {
                    setExisted(first, column);
                } else {
                    int lengthDiff = Math.abs(data[0].length - column.length);
                    if (lengthDiff > 0) {
                        insertNull(column, lengthDiff);
                    }
                    data[length++] = column;
                }
            }
        }
    }

    /*
     * Needed in a case we have nested list with maps - if so,
     * we get list of columns of different length
     */
    private void processListMap(List<Column> columns) {
        boolean nullsOnly = columns.stream().noneMatch(Objects::nonNull);
        if (nullsOnly) return;
        int minLength = columns.get(0).length;
        int maxLength = minLength;
        for (int i = 1; i < columns.size(); i++) {
            int length = columns.get(i).length;
            minLength = Math.min(minLength, length);
            maxLength = Math.max(maxLength, length);
        }
        if (maxLength - minLength == 0) return;
        for (Column column: columns) {
            appendNull(column, maxLength - column.length);
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
    public void set(int index, List<Column> value) {
        throw new UnsupportedOperationException("Not supported");
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
        if (Arrays.stream(data).allMatch(Objects::isNull)) {
            return returnArray;
        }
        System.arraycopy(data, 0, returnArray, 0, length);
        return returnArray;
    }

    /**
     * complexTypeConverter return always StringColumn for null values because it's not possible to detect type.
     * If not null value appears and its corresponding column is not StringColumn
     * and all previous values in old column are nulls - replace with new column type and insert nulls before value
     */
    private void setExisted(Column existed, Column newColumn) {
        boolean allNulls = true;
        for (int i = 0; i < existed.length; i++) {
            if (existed.get(i) != null) {
                allNulls = false;
                break;
            }
        }
        if (allNulls && !existed.getType().equals(newColumn.getType())) {
            insertNull(newColumn, existed.length - 1);
            for (int i = 0; i < data.length; i++) {
                if (data[i] == existed) {
                    data[i] = newColumn;
                    return;
                }
            }
        }
        for (int i = 0; i < newColumn.length; i++) {
            Object value = newColumn.get(i);
            if (existed.getType().equals(Types.BOOL)) {
                value = ((int) value) == 1;
            }
            if (i == 0) {
                existed.set(existed.length - 1, value); // replace added null
            } else {
                existed.add(value);
            }
        }
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            Column[] newData = new Column[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    private void insertNull(Column column, int nullCount) {
        if (column == null) return;
        try {
            Field data ;
            if (column.getType().equals(Types.BIG_INT)) {
                data = column.getClass().getSuperclass().getDeclaredField(DATA_FIELD_NAME);
            } else {
                data = column.getClass().getDeclaredField(DATA_FIELD_NAME);
            }
            data.setAccessible(true);
            Class<?> componentType = data.get(column).getClass().getComponentType();
            Object array = Array.newInstance(componentType.isPrimitive()
                    ? getWrapperClassForPrimitive(componentType) : componentType, column.length + nullCount);
            for (int j = nullCount, i = 0; i < column.length; j++, i++) {
                if (column.getType().equals(Types.BOOL)) {
                    boolean value = ((Integer) column.get(i)) == 1;
                    Array.set(array, j, value);
                } else {
                    Array.set(array, j, column.get(i));
                }
            }
            data.setAccessible(false);
            column.empty();
            column.addAll((convertToObjectArray(array)));
        } catch (NoSuchFieldException | IllegalAccessException e) {
            throw new RuntimeException("Something went wrong when inserting null", e);
        }
    }

    private Object[] convertToObjectArray(Object array) {
        Class ofArray = array.getClass().getComponentType();
        if (ofArray.isPrimitive()) {
            List ar = new ArrayList();
            int length = Array.getLength(array);
            for (int i = 0; i < length; i++) {
                ar.add(Array.get(array, i));
            }
            return ar.toArray();
        }
        else {
            return (Object[]) array;
        }
    }

    private Class<?> getWrapperClassForPrimitive(Class<?> primitiveClass) {
        switch (primitiveClass.getName()) {
            case "double":
                return Double.class;
            case "float":
                return Float.class;
            case "int":
                return Integer.class;
            default:
                throw new RuntimeException("Couldn't find wrapper");
        }
    }

    private void insertNull(Column[] columns, int nullCount) {
        for (Column column: columns) {
            insertNull(column, nullCount);
        }
    }

    private void appendNull(Column column, int nullCount) {
        if (column == null) return;
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
