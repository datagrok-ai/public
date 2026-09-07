package serialization;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.ZoneOffset;

public class DateTimeColumn extends AbstractColumn<Double> {
    private static final String TYPE = Types.DATE_TIME;
    private static final long MICROS_PER_DAY = 86_400_000_000L;
    private static final int MICROSECOND = 7;

    // dateTime:int (id 2) when its components beat dateTime:microseconds; false restores id 3 only.
    public static boolean COMPONENT_ENCODER =
            Boolean.parseBoolean(System.getProperty("grok.connect.dateTimeComponentEncoder", "true"));

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
        if (!COMPONENT_ENCODER || !encodeComponents(buf, false)) {
            buf.writeInt32(3);
            buf.writeFloat64List(data, 0, length);
        }
    }

    // Layout of DateTimeIntEncoder.encode (date_time_column_encoders.dart:228-258): per component a present flag,
    // then the widened values as a nested IntColumn payload; a component is absent only when all its values are 0
    // (_setDataValue allocates on the first non-zero), year/month/day are always present. Unless [force], nothing is
    // written and false is returned when the payload is not smaller than length * 8 or when a microsecond component
    // exists - dart2js rounds microseconds to milliseconds in DateTime.utc (core_patch.dart:256), so the browser
    // reader could not restore them.
    boolean encodeComponents(BufferAccessor buf, boolean force) {
        int[][] components = components();
        IntColumn.Encoding[] encodings = new IntColumn.Encoding[8];
        long payloadBytes = 8;
        for (int c = 0; c < 8; c++) {
            if (c >= 3 && allZero(components[c]))
                continue;
            encodings[c] = IntColumn.Encoding.choose(components[c], length);
            payloadBytes += 4 + encodings[c].size;
        }
        if (!force && (encodings[MICROSECOND] != null || payloadBytes >= length * 8L))
            return false;

        buf.writeInt32(2);
        for (int c = 0; c < 8; c++) {
            buf.writeInt8((byte) (encodings[c] != null ? 1 : 0));
            if (encodings[c] != null)
                encodings[c].write(buf, components[c], length);
        }
        return true;
    }

    // year, month, day, hour, minute, second, millisecond, microsecond per row, computed in UTC from the epoch-us
    // value like Dart's _getDataValue (toUtc) + _setDataValue. None rows carry (1,1,1,0,0,0,0,0), which every
    // reader maps back to None - as it does a genuine 0001-01-01T00:00:00Z (same ambiguity as the Dart writer).
    private int[][] components() {
        int[][] components = new int[8][length];
        for (int n = 0; n < length; n++) {
            if (data[n] == FloatColumn.None) {
                components[0][n] = 1;
                components[1][n] = 1;
                components[2][n] = 1;
                continue;
            }
            long micros = (long) data[n];
            LocalDate date = LocalDate.ofEpochDay(Math.floorDiv(micros, MICROS_PER_DAY));
            long ofDay = Math.floorMod(micros, MICROS_PER_DAY);
            components[0][n] = date.getYear();
            components[1][n] = date.getMonthValue();
            components[2][n] = date.getDayOfMonth();
            components[3][n] = (int) (ofDay / 3_600_000_000L);
            components[4][n] = (int) (ofDay / 60_000_000L % 60);
            components[5][n] = (int) (ofDay / 1_000_000L % 60);
            components[6][n] = (int) (ofDay / 1_000L % 1000);
            components[MICROSECOND][n] = (int) (ofDay % 1000);
        }
        return components;
    }

    private static boolean allZero(int[] values) {
        for (int v : values)
            if (v != 0)
                return false;
        return true;
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
