package serialization;

import java.util.*;
import java.lang.Math;
import java.nio.charset.*;


// A convenient class for binary data serialization.
@SuppressWarnings("unused")
public class BufferAccessor {
    // Internal data storage.
    private byte[] buf;

    // A view to the storage (convenient for casting content to different types).
    private ByteData view;

    // Current position in the buffer. All read/write methods use it.
    public int bufPos = 0;

    // Pending header to prepend on toUint8List().
    private byte[] pendingHeader;

    // Version of format.
    public static final String VERSION = "0.1.0-21a23d8e";

    public static final int TYPE = 0xA000;
    public static final int LIST = 0x0100;

    public static final int STRING = 1;
    public static final int FLOAT_32 = 2;
    public static final int UINT_8 = 3;
    public static final int UINT_16 = 4;
    public static final int UINT_32 = 5;
    public static final int INT_8 = 6;
    public static final int INT_16 = 7;
    public static final int INT_32 = 8;
    public static final int BOOL = 9;
    public static final int FLOAT_64 = 10;

    public static final int TYPE_COLUMN = TYPE + 40;
    public static final int TYPE_DATA_FRAME = TYPE + 41;
    public static final int TYPE_STRING_MAP = TYPE + 42;

    public static final int TYPE_FLOAT_32_LIST = TYPE + LIST + FLOAT_32;
    public static final int TYPE_FLOAT_64_LIST = TYPE + LIST + FLOAT_64;
    public static final int TYPE_UINT_8_LIST = TYPE + LIST + UINT_8;
    public static final int TYPE_UINT_16_LIST = TYPE + LIST + UINT_16;
    public static final int TYPE_UINT_32_LIST = TYPE + LIST + UINT_32;
    public static final int TYPE_INT_8_LIST = TYPE + LIST + INT_8;
    public static final int TYPE_INT_16_LIST = TYPE + LIST + INT_16;
    public static final int TYPE_INT_32_LIST = TYPE + LIST + INT_32;
    public static final int TYPE_STRING_LIST = TYPE + LIST + STRING;
    // Writes two bytes that determine the type of the entity that follows.
    private void writeTypeCode(int typeCode) {
        writeInt16((short) typeCode);
    }

    public BufferAccessor() {
        buf = new byte[100];
        view = new ByteData(buf);
    }

    public BufferAccessor(byte[] list) {
        buf = list;
        view = new ByteData(buf);
    }

    void writeString(String value) {
        byte[] bytes = value == null ? null : value.getBytes(StandardCharsets.UTF_8);
        writeUint8List(bytes);
    }

    void writeInt8(byte value) {
        _ensureSpace(1);
        view.setInt8(bufPos, value);
        bufPos += 1;
    }

    void writeInt16(short value) {
        _ensureSpace(2);
        view.setInt16(bufPos, value);
        bufPos += 2;
    }

    public void writeInt32(int value) {
        _ensureSpace(4);
        view.setInt32(bufPos, value);
        bufPos += 4;
    }

    public void writeInt64(long value) {
        _ensureSpace(8);
        view.setFloat64(bufPos, (double) value);
        bufPos += 8;
    }

    public void writeFloat32List(float[] values, int... idxs) {
        int start = (idxs.length > 0) ? idxs[0] : 0;
        int count = (idxs.length > 1) ? idxs[1] : values.length - start;
        beginList(TYPE_FLOAT_32_LIST, 4, count);
        for (int i = 0; i < count; i++)
            view.setFloat32(bufPos + i * 4, values[start + i]);
        bufPos += count * 4;
    }

    public void writeFloat64List(double[] values, int... idxs) {
        int start = (idxs.length > 0) ? idxs[0] : 0;
        int count = (idxs.length > 1) ? idxs[1] : values.length - start;
        beginList(TYPE_FLOAT_64_LIST, 8, count);
        for (int i = 0; i < count; i++)
            view.setFloat64(bufPos + i * 8, values[start + i]);
        bufPos += count * 8;
    }

    void writeUint8List(byte[] values, int... idxs) {
        writeTypeCode(TYPE_UINT_8_LIST);
        if (values == null) {
            writeInt64(-1);
            return;
        }

        int start = (idxs.length > 0) ? idxs[0] : 0;
        int count = (idxs.length > 1) ? idxs[1] : values.length - start;

        _ensureSpace(8 + count);

        writeInt64(count);
        System.arraycopy(values, start, buf, bufPos, count);
        bufPos += count;
    }

    public void writeInt32List(int[] values, int... idxs) {
        writeInt32Array(values, TYPE_INT_32_LIST, idxs);
    }

    public void writeUint32List(int[] values, int... idxs) {
        writeInt32Array(values, TYPE_UINT_32_LIST, idxs);
    }

    private void writeInt32Array(int[] values, int typeCode, int... idxs) {
        int start = (idxs.length > 0) ? idxs[0] : 0;
        int count = (idxs.length > 1) ? idxs[1] : values.length - start;
        beginList(typeCode, 4, count);
        for (int i = 0; i < count; i++)
            view.setInt32(bufPos + i * 4, values[start + i]);
        bufPos += count * 4;
    }

    private void beginList(int typeCode, int elementSize, int count) {
        writeTypeCode(typeCode);
        _ensureSpace(8 + count * elementSize);
        writeInt64(count);
    }

    public void writeStringList(String[] values, int... idxs) {
        writeTypeCode(TYPE_STRING_LIST);
        int start = (idxs.length > 0) ? idxs[0] : 0;
        int count = (idxs.length > 1) ? idxs[1] : values.length - start;
        writeTypeCode(TYPE_UINT_8_LIST);
        int begin = bufPos;
        int[] lengths = new int[count];
        writeInt64(0);
        for (int n = 0; n < count; n++) {
            String str = values[start + n];
            if (str != null) {
                lengths[n] = str.length();
                byte[] bytes = str.getBytes(StandardCharsets.UTF_8);
                _ensureSpace(bytes.length);
                System.arraycopy(bytes, 0, buf, bufPos, bytes.length);
                bufPos += bytes.length;
            } else
                lengths[n] = -1;
        }
        int end = bufPos;
        bufPos = begin;
        writeInt64(end - begin - 8);
        bufPos = end;
        writeInt32List(lengths);
    }

    void writeColumn(Column<?> col) {
        if (col instanceof ComplexTypeColumn) {
            col.encode(this);
            return;
        }
        writeTypeCode(TYPE_COLUMN);
        writeString(col.getName());
        writeString(col.getType());
        writeStringMap(col.getTags());
        col.encode(this);
    }

    // Serializes a [DataFrame]. It can be null. Returns columns offsets [offsets].
    public int[] writeDataFrame(DataFrame dataFrame) {
        writeTypeCode(TYPE_DATA_FRAME);
        writeInt64((dataFrame.rowCount == null) ? -1 : dataFrame.rowCount);

        int flatCount = 0;
        for (int n = 0; n < dataFrame.getColumnCount(); n++) {
            Column<?> col = dataFrame.getColumn(n);
            flatCount += (col instanceof ComplexTypeColumn)
                ? ((ComplexTypeColumn) col).getEncodedColumnCount()
                : 1;
        }
        writeInt64(flatCount);

        writeString(dataFrame.name);
        writeStringMap(dataFrame.getTags());
        int[] offsets = new int[dataFrame.getColumnCount()];
        for (int n = 0; n < dataFrame.getColumnCount(); n++) {
            offsets[n] = bufPos;
            writeColumn(dataFrame.getColumn(n));
        }

        return offsets;
    }

    private void writeStringMap(Map<String, String> map) {
        writeTypeCode(TYPE_STRING_MAP);

        if (map == null) {
            writeInt32(-1);
            return;
        }

        writeInt32(map.size());
        for (String key : map.keySet()) {
            writeString(key);
            writeString(map.get(key));
        }
    }

    private void _ensureSpace(int extraLength) {
        if (bufPos + extraLength > buf.length) {
            byte[] newBuf = new byte[buf.length * 2 + Math.max(0, bufPos + extraLength - buf.length * 2)];
            System.arraycopy(buf, 0, newBuf, 0, buf.length);
            buf = newBuf;
            view = new ByteData(buf);
        }
    }

    public byte[] toUint8List() {
        if (pendingHeader != null) {
            byte[] result = new byte[pendingHeader.length + bufPos];
            System.arraycopy(pendingHeader, 0, result, 0, pendingHeader.length);
            System.arraycopy(buf, 0, result, pendingHeader.length, bufPos);
            pendingHeader = null;
            return result;
        }
        byte[] result = new byte[bufPos];
        System.arraycopy(buf, 0, result, 0, bufPos);
        return result;
    }

    // Prepends a string header. The header is stored separately and
    // concatenated with the main buffer on the next toUint8List() call,
    // avoiding an O(n) shift of the entire buffer.
    public void insertStringHeader(String header) {
        BufferAccessor headerBuf = new BufferAccessor();
        headerBuf.writeString(header);
        pendingHeader = headerBuf.toUint8List();
    }

    // Writes raw bytes to buffer with no overhead.
    void writeRawBytes(byte[] bytes) {
        _ensureSpace(bytes.length);
        System.arraycopy(bytes, 0, buf, bufPos, bytes.length);
        bufPos += bytes.length;
    }

    // ===================================================================
    // Read half (mirror of the write methods; ports serialization.dart).
    // ===================================================================

    private static void _checkTypeCode(int typeCode) {
        if ((typeCode & 0xF000) != TYPE)
            throw new RuntimeException("Invalid type identifier: " + typeCode);
    }

    // Reads the previously saved type code and asserts its equality to [expectedTypeCode].
    public void readTypeCode(int expectedTypeCode) {
        int code = readUint16();
        _checkTypeCode(code);
        if (code != expectedTypeCode)
            throw new RuntimeException("Type code mismatch. Read code " + code
                    + ". Expected code " + expectedTypeCode);
    }

    public int readInt8() {
        return view.getInt8(bufPos++);
    }

    public int readInt16() {
        short v = view.getInt16(bufPos);
        bufPos += 2;
        return v;
    }

    public int readInt32() {
        int v = view.getInt32(bufPos);
        bufPos += 4;
        return v;
    }

    public int readUint8() {
        return view.getUint8(bufPos++);
    }

    public int readUint16() {
        int v = view.getUint16(bufPos);
        bufPos += 2;
        return v;
    }

    // NOTE dart2js has no Int64 accessor, so int64 is encoded as a float64.
    // Reader reads a double, rejects values above 2^53, and casts.
    public long readInt64() {
        double f = view.getFloat64(bufPos);
        bufPos += 8;
        if (f > 0x1FFFFFFFFFFFFFL)
            throw new RuntimeException("Value too big");
        return (long) f;
    }

    public float readFloat32() {
        float v = view.getFloat32(bufPos);
        bufPos += 4;
        return v;
    }

    public double readFloat64() {
        double v = view.getFloat64(bufPos);
        bufPos += 8;
        return v;
    }

    public String readString() {
        byte[] bytes = readUint8List();
        return bytes == null ? null : new String(bytes, StandardCharsets.UTF_8);
    }

    public byte[] readUint8List() {
        readTypeCode(TYPE_UINT_8_LIST);
        int len = (int) readInt64();
        if (len == -1)
            return null;
        byte[] list = new byte[len];
        System.arraycopy(buf, bufPos, list, 0, len);
        bufPos += len;
        return list;
    }

    public byte[] readInt8List() {
        readTypeCode(TYPE_INT_8_LIST);
        int len = (int) readInt64();
        if (len == -1)
            return null;
        byte[] list = new byte[len];
        System.arraycopy(buf, bufPos, list, 0, len);
        bufPos += len;
        return list;
    }

    public short[] readInt16List() {
        readTypeCode(TYPE_INT_16_LIST);
        int len = (int) readInt64();
        short[] list = new short[len];
        for (int i = 0; i < len; i++)
            list[i] = view.getInt16(bufPos + i * 2);
        bufPos += len * 2;
        return list;
    }

    public short[] readUint16List() {
        readTypeCode(TYPE_UINT_16_LIST);
        int len = (int) readInt64();
        short[] list = new short[len];
        for (int i = 0; i < len; i++)
            list[i] = view.getInt16(bufPos + i * 2);
        bufPos += len * 2;
        return list;
    }

    public int[] readInt32List() {
        readTypeCode(TYPE_INT_32_LIST);
        int len = (int) readInt64();
        int[] list = new int[len];
        for (int i = 0; i < len; i++)
            list[i] = view.getInt32(bufPos + i * 4);
        bufPos += len * 4;
        return list;
    }

    public int[] readUint32List() {
        readTypeCode(TYPE_UINT_32_LIST);
        int len = (int) readInt64();
        int[] list = new int[len];
        for (int i = 0; i < len; i++)
            list[i] = view.getInt32(bufPos + i * 4);
        bufPos += len * 4;
        return list;
    }

    public float[] readFloat32List() {
        readTypeCode(TYPE_FLOAT_32_LIST);
        int len = (int) readInt64();
        float[] list = new float[len];
        for (int i = 0; i < len; i++)
            list[i] = view.getFloat32(bufPos + i * 4);
        bufPos += len * 4;
        return list;
    }

    public double[] readFloat64List() {
        readTypeCode(TYPE_FLOAT_64_LIST);
        int len = (int) readInt64();
        double[] list = new double[len];
        for (int i = 0; i < len; i++)
            list[i] = view.getFloat64(bufPos + i * 8);
        bufPos += len * 8;
        return list;
    }

    // String list. Lengths are UTF-16 code-unit counts (NOT UTF-8 byte lengths):
    // decode the whole byte region as one UTF-8 string, then slice by cumulative
    // code-unit offsets (String.substring is UTF-16 in Java, as in Dart).
    // A -1 length decodes to '' (also the string-column None slot).
    public String[] readStringList() {
        readTypeCode(TYPE_STRING_LIST);
        String str = readString();
        int[] lengths = readInt32List();
        String[] list = new String[lengths.length];
        int start = 0;
        for (int n = 0; n < lengths.length; n++) {
            int len = lengths[n];
            if (len >= 0) {
                list[n] = str.substring(start, start + len);
                start += len;
            } else
                list[n] = "";
        }
        return list;
    }

    public Map<String, String> readStringMap() {
        readTypeCode(TYPE_STRING_MAP);
        int pairs = readInt32();
        if (pairs == -1)
            return null;
        Map<String, String> map = new LinkedHashMap<>();
        for (int i = 0; i < pairs; i++) {
            String key = readString();
            map.put(key, readString());
        }
        return map;
    }

    public byte[] readRawBytes(int length) {
        byte[] list = new byte[length];
        System.arraycopy(buf, bufPos, list, 0, length);
        bufPos += length;
        return list;
    }

    // Reads a Column: TYPE_COLUMN, name, typeName, tags, then encoder payload.
    public Column<?> readColumn() {
        readTypeCode(TYPE_COLUMN);
        String colName = readString();
        String colType = readString();
        Map<String, String> tags = readStringMap();

        Column<?> col = createColumn(colType, colName);
        col.decode(this);
        if (tags != null && !tags.isEmpty())
            col.getTags().putAll(tags);

        return col;
    }

    private static Column<?> createColumn(String type, String name) {
        switch (type) {
            case Types.INT:       return new IntColumn(name, 0);
            case Types.FLOAT:     return new FloatColumn(name, 0);
            case Types.BOOL:      return new BoolColumn(name, 0);
            case Types.STRING:    return new StringColumn(name, 0);
            case Types.DATE_TIME: return new DateTimeColumn(name, 0);
            case Types.BIG_INT:   return new BigIntColumn(name, 0);
            default:
                throw new RuntimeException("Unknown column type for decoding: " + type);
        }
    }

    // Reads a DataFrame, which can be null. When [offsets] is provided, seeks to
    // each column's absolute offset before reading it.
    public DataFrame readDataFrame(int[] offsets) {
        readTypeCode(TYPE_DATA_FRAME);
        long rowCount = readInt64();
        if (rowCount == -1)
            return null;

        int colCount = (int) readInt64();
        String name = readString();
        Map<String, String> tags = readStringMap();

        DataFrame df = new DataFrame();
        df.name = name;
        if (tags != null && !tags.isEmpty())
            df.getTags().putAll(tags);

        for (int i = 0; i < colCount; i++) {
            if (offsets != null)
                bufPos = offsets[i];
            df.addColumn(readColumn());
        }
        df.rowCount = (int) rowCount;
        return df;
    }
}
