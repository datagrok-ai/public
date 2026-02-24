package serialization;

public interface Column<T> extends Taggable {
    String getName();
    void setName(String name);
    int getLength();
    String getType();
    void encode(BufferAccessor buf);
    void add(T value);
    void addAll(T[] value);
    T get(int idx);
    void set(int index, T value);
    long memoryInBytes();
    boolean isNone(int idx);
    void empty();
    Object toArray();
    int getInitColumnSize();
    void setInitColumnSize(int size);

    static Column<?> getColumnForType(String type, String name, int initColumnSize) {
        switch (type) {
            case Types.FLOAT:
                return new FloatColumn(name, initColumnSize);
            case Types.INT:
                return new IntColumn(name, initColumnSize);
            case Types.BIG_INT:
                return new BigIntColumn(name, initColumnSize);
            case Types.DATE_TIME:
                return new DateTimeColumn(name, initColumnSize);
            case Types.COLUMN_LIST:
                return new ComplexTypeColumn(name, initColumnSize);
            case Types.BOOL:
                return new BoolColumn(name, initColumnSize);
            case Types.STRING:
            default:
                return new StringColumn(name, initColumnSize);
        }
    }
}
