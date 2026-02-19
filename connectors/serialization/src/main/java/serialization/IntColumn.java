package serialization;

public class IntColumn extends AbstractColumn<Integer> {
    private static final String TYPE = Types.INT;
    public static final int None = -2147483648;

    private int[] data;

    public IntColumn(String name) {
        super(name);
        data = new int[initColumnSize];
    }

    public IntColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        data = new int[initColumnSize];
    }

    public IntColumn(String name, Integer[] values) {
        super(name);
        data = new int[initColumnSize];
        addAll(values);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        data = new int[initColumnSize];
    }

    @Override
    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);
        buf.writeInt8((byte)0);
        buf.writeInt32List(data, 0, length);
    }

    @Override
    public void add(Integer value) {
        ensureSpace(1);
        data[length++] = (value != null) ? value : None;
    }

    @Override
    public void addAll(Integer[] values) {
        ensureSpace(values.length);
        for (Integer value : values)
            data[length++] = (value != null) ? value : None;
    }

    @Override
    public Integer get(int idx) {
        return data[idx];
    }

    @Override
    public void set(int index, Integer value) {
        data[index] = value;
    }

    @Override
    public long memoryInBytes() {
        return (long) data.length * 4;
    }

    @Override
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

    @Override
    public Object toArray() {
        return data;
    }
}
