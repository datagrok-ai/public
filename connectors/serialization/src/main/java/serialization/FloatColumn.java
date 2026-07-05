package serialization;

public class FloatColumn extends AbstractColumn<Float> {
    private static final String TYPE = Types.FLOAT;
    static final double None = 2.6789344063684636e-34;

    // Single-precision storage (float:raw, encoder id 1).
    private float[] data;

    // Double-precision storage (float:raw64, encoder id 5). Dart's default float
    // encoder writes float64, so Datlas-written double columns arrive here.
    private double[] data64;
    private boolean doublePrecision;

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

    // Builds a double-precision float column (encoder id 5).
    public static FloatColumn double64(String name, Double[] values) {
        FloatColumn c = new FloatColumn(name, Math.max(values == null ? 0 : values.length, 1));
        c.doublePrecision = true;
        c.data = null;
        c.data64 = new double[c.initColumnSize];
        c.length = 0;
        if (values != null)
            for (Double value : values)
                c.addDouble(value);
        return c;
    }

    public boolean isDoublePrecision() {
        return doublePrecision;
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        if (doublePrecision)
            data64 = new double[initColumnSize];
        else
            data = new float[initColumnSize];
    }

    @Override
    public void encode(BufferAccessor buf) {
        if (doublePrecision) {
            buf.writeInt32(5);       // Encoder ID (float:raw64)
            buf.writeInt16((short)1); // metadata param count (excluding archive)
            buf.writeInt8((byte)1);   // doublePrecision
            buf.writeInt8((byte)0);   // Archive
            buf.writeFloat64List(data64, 0, length);
        } else {
            buf.writeInt32(1);       // Encoder ID (float:raw)
            buf.writeInt8((byte)0);   // Archive
            buf.writeFloat32List(data, 0, length);
        }
    }

    @Override
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        switch (id) {
            case 1: // float:raw (legacy, no metadata header) - float32
                doublePrecision = false;
                if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
                    data = ByteData.toFloat32List(Zlib.inflate(buf.readUint8List()));
                else
                    data = buf.readFloat32List();
                length = data.length;
                break;
            case 2: // float:rle (legacy, no metadata header) - float32
                doublePrecision = false;
                data = serialization.codecs.FloatRle.decode32(buf);
                length = data.length;
                break;
            case 3: // float:fcp - float32 only
                doublePrecision = false;
                data = new serialization.codecs.FloatFcp().decode(buf);
                length = data.length;
                break;
            case 4: // float:rle64 (metadata header)
                buf.readInt16();
                doublePrecision = buf.readInt8() == 1;
                if (doublePrecision) {
                    data64 = serialization.codecs.FloatRle.decode64(buf);
                    length = data64.length;
                } else {
                    data = serialization.codecs.FloatRle.decode32(buf);
                    length = data.length;
                }
                break;
            case 5: // float:raw64 (metadata header)
                buf.readInt16();
                doublePrecision = buf.readInt8() == 1;
                if (doublePrecision) {
                    if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
                        data64 = ByteData.toFloat64List(Zlib.inflate(buf.readUint8List()));
                    else
                        data64 = buf.readFloat64List();
                    length = data64.length;
                } else {
                    if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
                        data = ByteData.toFloat32List(Zlib.inflate(buf.readUint8List()));
                    else
                        data = buf.readFloat32List();
                    length = data.length;
                }
                break;
            default:
                throw new RuntimeException("decoding " + name + ": float encoder " + id + " not found");
        }
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

    public void addDouble(Double value) {
        ensureSpace64(1);
        data64[length++] = (value != null) ? value : None;
    }

    @Override
    public Float get(int idx) {
        return doublePrecision ? (float) data64[idx] : data[idx];
    }

    // Returns the value at full stored precision.
    public double getDouble(int idx) {
        return doublePrecision ? data64[idx] : data[idx];
    }

    @Override
    public void set(int index, Float value) {
        if (doublePrecision)
            data64[index] = value;
        else
            data[index] = value;
    }

    @Override
    public long memoryInBytes() {
        return doublePrecision ? (long) data64.length * 8 : (long) data.length * 4;
    }

    @Override
    public boolean isNone(int idx) {
        return doublePrecision ? data64[idx] == None : data[idx] == None;
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            float[] newData = new float[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    private void ensureSpace64(int extraLength) {
        if (length + extraLength > data64.length) {
            double[] newData = new double[data64.length * 2 + Math.max(0, length + extraLength - data64.length * 2)];
            System.arraycopy(data64, 0, newData, 0, data64.length);
            data64 = newData;
        }
    }

    @Override
    public Object toArray() {
        return doublePrecision ? (Object) data64 : (Object) data;
    }
}
