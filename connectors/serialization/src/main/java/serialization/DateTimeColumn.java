package serialization;

import java.time.LocalDateTime;
import java.time.ZoneOffset;

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
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        switch (id) {
            case 1: // dateTime:raw - 8 optional component arrays.
                decodeComponents(buf, false);
                break;
            case 2: // dateTime:int - components as nested int columns.
                decodeComponents(buf, true);
                break;
            case 3: // dateTime:microseconds - epoch-us float64. None = FloatColumn.None.
                data = buf.readFloat64List();
                length = data.length;
                break;
            default:
                throw new RuntimeException("decoding " + name + ": datetime encoder " + id + " not found");
        }
    }

    // Ports date_time_column_encoders.dart:23-64 (raw) and :186-221 (int):
    // reads the 8 component arrays and materializes each row to epoch-us like
    // Dart :61-63 (dtu.create -> DateTime.utc -> microsecondsSinceEpoch).
    private void decodeComponents(BufferAccessor buf, boolean asIntColumns) {
        short[] year = asIntColumns ? readInt16IntCol(buf) : readInt16Component(buf);
        byte[] month = asIntColumns ? readInt8IntCol(buf) : readInt8Component(buf);
        byte[] day = asIntColumns ? readInt8IntCol(buf) : readInt8Component(buf);
        byte[] hour = asIntColumns ? readInt8IntCol(buf) : readInt8Component(buf);
        byte[] minute = asIntColumns ? readInt8IntCol(buf) : readInt8Component(buf);
        byte[] second = asIntColumns ? readInt8IntCol(buf) : readInt8Component(buf);
        short[] millisecond = asIntColumns ? readInt16IntCol(buf) : readInt16Component(buf);
        short[] microsecond = asIntColumns ? readInt16IntCol(buf) : readInt16Component(buf);

        length = year.length;
        data = new double[length];
        for (int n = 0; n < length; n++) {
            if (isNoneComponent(year, month, day, hour, minute, second, millisecond, microsecond, n))
                data[n] = FloatColumn.None;
            else
                data[n] = componentToMicros(year, month, day, hour, minute, second, millisecond, microsecond, n);
        }
    }

    // Ports the raw-encoder readInt16List/readInt8List closures (:24-50).
    private static short[] readInt16Component(BufferAccessor buf) {
        if (buf.readInt8() != 1)
            return null;
        if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
            return ByteData.toInt16List(Zlib.inflate(buf.readUint8List()));
        return buf.readInt16List();
    }

    private static byte[] readInt8Component(BufferAccessor buf) {
        if (buf.readInt8() != 1)
            return null;
        if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
            return Zlib.inflate(buf.readUint8List());
        return buf.readInt8List();
    }

    // Ports the int-encoder decodeInt16List/decodeInt8List closures (:187-207):
    // a present flag then a nested int column (encoder decode, no TYPE_COLUMN header).
    private static short[] readInt16IntCol(BufferAccessor buf) {
        if (buf.readInt8() == 0)
            return null;
        IntColumn ic = new IntColumn("", 0);
        ic.decode(buf);
        int[] d = (int[]) ic.toArray();
        short[] r = new short[d.length];
        for (int i = 0; i < d.length; i++)
            r[i] = (short) d[i];
        return r;
    }

    private static byte[] readInt8IntCol(BufferAccessor buf) {
        if (buf.readInt8() == 0)
            return null;
        IntColumn ic = new IntColumn("", 0);
        ic.decode(buf);
        int[] d = (int[]) ic.toArray();
        byte[] r = new byte[d.length];
        for (int i = 0; i < d.length; i++)
            r[i] = (byte) d[i];
        return r;
    }

    // Ports _isNone (date_time_column_encoders.dart:106-112).
    private static boolean isNoneComponent(short[] year, byte[] month, byte[] day,
            byte[] hour, byte[] minute, byte[] second, short[] millisecond, short[] microsecond, int idx) {
        return year[idx] == 1 && month[idx] == 1 && day[idx] == 1
                && (hour == null || hour[idx] == 0)
                && (minute == null || minute[idx] == 0)
                && (second == null || second[idx] == 0)
                && (millisecond == null || millisecond[idx] == 0)
                && (microsecond == null || microsecond[idx] == 0);
    }

    // Ports _getDataValue -> microsecondsSinceEpoch (dtu.create = DateTime.utc).
    private static double componentToMicros(short[] year, byte[] month, byte[] day,
            byte[] hour, byte[] minute, byte[] second, short[] millisecond, short[] microsecond, int idx) {
        int ms = millisecond != null ? millisecond[idx] : 0;
        int us = microsecond != null ? microsecond[idx] : 0;
        int nanoOfSecond = (ms * 1000 + us) * 1000;
        LocalDateTime ldt = LocalDateTime.of(
                year[idx],
                month[idx],
                day[idx],
                hour != null ? hour[idx] : 0,
                minute != null ? minute[idx] : 0,
                second != null ? second[idx] : 0,
                nanoOfSecond);
        long micros = ldt.toEpochSecond(ZoneOffset.UTC) * 1_000_000L + ldt.getNano() / 1000L;
        return (double) micros;
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
