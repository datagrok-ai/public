package serialization;

import java.util.ArrayList;
import java.util.List;

// Read-capable minimal port of ByteArrayColumn (byteArray:raw, encoder id 1).
// Storage is List<byte[]>; a null entry is None. Ports
// byte_array_column_encoders.dart:11-22. The Java writer never produced these and
// the mutation transport does not use them - the reader must not choke on them.
public class ByteArrayColumn extends AbstractColumn<byte[]> {
    private static final String TYPE = Types.BYTE_ARRAY;

    private List<byte[]> data = new ArrayList<>();

    public ByteArrayColumn(String name) {
        super(name);
    }

    public ByteArrayColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        data = new ArrayList<>();
    }

    @Override
    public void encode(BufferAccessor buf) {
        throw new UnsupportedOperationException("ByteArrayColumn encode is not supported");
    }

    @Override
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        if (id != 1)
            throw new RuntimeException("decoding " + name + ": byteArray encoder " + id + " not found");
        if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
            deserialize(new BufferAccessor(Zlib.inflate(buf.readUint8List())));
        else
            deserialize(buf);
    }

    private void deserialize(BufferAccessor buf) {
        int colLength = (int) buf.readInt64();
        data = new ArrayList<>(Math.max(colLength, 0));
        for (int i = 0; i < colLength; i++)
            data.add(buf.readUint8List());
        length = colLength;
    }

    @Override
    public void add(byte[] value) {
        data.add(value);
        length++;
    }

    @Override
    public void addAll(byte[][] values) {
        for (byte[] value : values)
            add(value);
    }

    @Override
    public byte[] get(int idx) {
        return data.get(idx);
    }

    @Override
    public void set(int index, byte[] value) {
        data.set(index, value);
    }

    @Override
    public long memoryInBytes() {
        long sum = 0;
        for (byte[] b : data)
            sum += (b == null) ? 0 : b.length;
        return sum;
    }

    @Override
    public boolean isNone(int idx) {
        return data.get(idx) == null;
    }

    @Override
    public Object toArray() {
        return data.toArray(new byte[0][]);
    }
}
