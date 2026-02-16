package serialization;

public class FloatColumn extends AbstractColumn<Float> {
    private static final String TYPE = Types.FLOAT;
    static final double None = 2.6789344063684636e-34;

    private float[] data;

    public FloatColumn(String name) {
        super(name);
        data = new float[initColumnSize];
    }

    public FloatColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        data = new float[initColumnSize];
    }

    public FloatColumn(String name, Float[] values) {
        super(name);
        data = new float[initColumnSize];
        addAll(values);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        data = new float[initColumnSize];
    }

    @Override
    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);  // Encoder ID
        buf.writeInt8((byte)0);   // Archive
        buf.writeFloat32List(data, 0, length);
    }

    @Override
    public void add(Float value) {
        ensureSpace(1);
        data[length++] = (value != null) ? value : (float)None;
    }

    @Override
    public void addAll(Float[] values) {
        ensureSpace(values.length);
        for (Float value : values)
            data[length++] = (value != null) ? value : (float)None;
    }

    @Override
    public Float get(int idx) {
        return data[idx];
    }

    @Override
    public void set(int index, Float value) {
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
            float[] newData = new float[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    @Override
    public Object toArray() {
        return data;
    }
}
