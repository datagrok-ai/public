package serialization;

import java.util.HashMap;
import java.util.Map;

// Column.
public abstract class Column<T> {
    public String name;
    public int length = 0;
    public Map<String, String> tags = new HashMap<>();
    public int initColumnSize = 100;

    public Column() {
    }

    public abstract String getType();
    public abstract void encode(BufferAccessor buf);
    public abstract void add(T value);
    public abstract void addAll(T[] value);
    public abstract Object get(int idx);
    public abstract void set(int index, T value);
    public abstract long memoryInBytes();
    public abstract boolean isNone(int idx);
    public abstract void empty();

    public static Column getColumnForType(String type, int initColumnSize) {
        switch (type) {
            case Types.FLOAT:
                return new FloatColumn(initColumnSize);
            case Types.INT:
                return new IntColumn(initColumnSize);
            case Types.BIG_INT:
                return new BigIntColumn(initColumnSize);
            case Types.DATE_TIME:
                return new DateTimeColumn(initColumnSize);
            case Types.COLUMN_LIST:
                return new ComplexTypeColumn();
            case Types.BOOL:
                return new BoolColumn(initColumnSize);
            case Types.STRING:
            default:
                return new StringColumn(initColumnSize);
        }
    }
}
