package serialization;

import java.util.*;


// Column.
public abstract class Column<T> {
    public String name;
    public int length = 0;
    public Map<String, String> tags = new HashMap<>();

    public Column() {
    }

    public abstract String getType();
    public abstract void encode(BufferAccessor buf);
    public abstract void add(T value);
    public abstract void addAll(T[] value);
    public abstract Object get(int idx);
    public abstract long memoryInBytes();
    public abstract boolean isNone(int idx);
}
