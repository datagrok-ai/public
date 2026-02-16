package serialization;

public class DateTimeColumn extends AbstractColumn<Double> {
    private static final String TYPE = Types.DATE_TIME;

    private double[] data;

    public DateTimeColumn(String name) {
        super(name);
        data = new double[initColumnSize];
    }

    public DateTimeColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        data = new double[initColumnSize];
    }

    public DateTimeColumn(String name, Double[] values) {
        super(name);
        data = new double[initColumnSize];
        addAll(values);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        data = new double[initColumnSize];
    }

    @Override
    public void encode(BufferAccessor buf) {
        buf.writeInt32(3); // Encoder ID
        buf.writeFloat64List(data, 0, length);
    }

    @Override
    public void add(Double value) {
        ensureSpace(1);
        setValue(length++, (value != null) ? value : FloatColumn.None);
    }

    @Override
    public void addAll(Double[] values) {
        ensureSpace(values.length);
        for (Double value : values)
            setValue(length++, (value != null) ? value : FloatColumn.None);
    }

    @Override
    public Double get(int idx) {
        return data[idx];
    }

    @Override
    public void set(int index, Double value) {
        data[index] = value;
    }

    @Override
    public long memoryInBytes() {
        return (long) data.length * 8;
    }

    @Override
    public boolean isNone(int idx) {
        return data[idx] == FloatColumn.None;
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            double[] newData = new double[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    private void setValue(int idx, Double value) {
        data[idx] = value;
    }

    @Override
    public Object toArray() {
        return data;
    }
}
