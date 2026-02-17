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

    private void writeInt16(short value) {
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
}
