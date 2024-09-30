package serialization;

// Integer column.
public class IntColumn extends Column<Integer> {
    private static final String TYPE = Types.INT;
    public static final int None = -2147483648;

    private int[] data;

    public IntColumn() {
        data = new int[initColumnSize];
    }

    public IntColumn(String name) {
        this();
        this.name = name;
    }

    public IntColumn(int initColumnSize) {
        this.initColumnSize = initColumnSize;
        data = new int[initColumnSize];
    }

    public IntColumn(Integer[] values) {
        data = new int[initColumnSize];
        addAll(values);
    }

    public String getType() {
        return TYPE;
    }

    public void empty() {
        length = 0;
        data = new int[initColumnSize];
    }

    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);
        buf.writeInt8((byte)0);
        buf.writeInt32List(data, 0, length);
    }

    public void add(Integer value) {
        ensureSpace(1);
        data[length++] = (value != null) ? value : None;
    }

    public void addAll(Integer[] values) {
        ensureSpace(values.length);
        for (int n = 0; n < values.length; n++)
            data[length++] = (values[n] != null) ? values[n] : None;
    }

    public Object get(int idx) {
        return data[idx];
    }

    @Override
    public void set(int index, Integer value) {
        data[index] = value;
    }

    @Override
    public long memoryInBytes() {
        return data.length * 4;
    }

    public boolean isNone(int idx) {
        return data[idx] == None;
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            int[] newData = new int[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    public int[] getData() {
        return data;
    }
}
