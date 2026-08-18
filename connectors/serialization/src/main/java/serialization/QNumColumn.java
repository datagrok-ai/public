package serialization;

// Read-capable minimal port of QNumColumn (qnum:raw, encoder id 1). Storage is a
// packed Float64List where the two least-significant mantissa bits hold the
// qualifier (<, =, >); see ddt/lib/src/utils/qnum.dart. The Java writer never
// produced qnum columns and the mutation transport does not use them - the reader
// only needs to not choke on them and expose values/none faithfully.
public class QNumColumn extends AbstractColumn<Double> {
    private static final String TYPE = Types.QNUM;
    static final double None = FloatColumn.None;

    private double[] data;

    public QNumColumn(String name) {
        super(name);
        data = new double[initColumnSize];
    }

    public QNumColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        data = new double[initColumnSize];
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
        throw new UnsupportedOperationException("QNumColumn encode is not supported");
    }

    @Override
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        if (id != 1)
            throw new RuntimeException("decoding " + name + ": qnum encoder " + id + " not found");
        if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
            data = ByteData.toFloat64List(Zlib.inflate(buf.readUint8List()));
        else
            data = buf.readFloat64List();
        length = data.length;
    }

    @Override
    public void add(Double value) {
        ensureSpace(1);
        data[length++] = (value != null) ? value : None;
    }

    @Override
    public void addAll(Double[] values) {
        ensureSpace(values.length);
        for (Double value : values)
            data[length++] = (value != null) ? value : None;
    }

    // Full stored value (with the qualifier bits still packed in).
    public double getRaw(int idx) {
        return data[idx];
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

    // Mirrors QNumColumn.isNone: strip the 2 qualifier bits, compare to None; NaN is none.
    @Override
    public boolean isNone(int idx) {
        double stripped = Double.longBitsToDouble(Double.doubleToRawLongBits(data[idx]) & ~0x3L);
        return stripped == None || Double.isNaN(data[idx]);
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            double[] newData = new double[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    @Override
    public Object toArray() {
        return data;
    }
}
