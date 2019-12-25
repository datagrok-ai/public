package serialization;


// Bool column.
public class BoolColumn extends Column<Boolean> {
    private static final String TYPE = Types.BOOL;

    private int[] data;

    public String getType() {
        return TYPE;
    }

    public BoolColumn() {
        data = new int[100];
    }

    public BoolColumn(Boolean[] values) {
        data = new int[100];
        addAll(values);
    }

    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);
        buf.writeInt64(length);
        buf.writeInt8((byte)0);
        buf.writeUint32List(data, 0, ((length + 0x1F) / 0x20));
    }

    public void add(Boolean value) {
        ensureSpace(1);
        if ((value != null) && value)
            data[length / 0x20] |= 1 << ((length % 0x20) & 0x1F);
        length++;
    }

    public void addAll(Boolean[] values) {
        ensureSpace(values.length);
        for (int n = 0; n < values.length; n++) {
            if ((values[n] != null) && values[n])
                data[length / 0x20] |= 1 << ((length % 0x20) & 0x1F);
            length++;
        }
    }

    @Override
    public long memoryInBytes() {
        return data.length * 4;
    }

    private void ensureSpace(int extraLength) {
        int lengthInInts = ((length + extraLength + 0x1F) / 0x20);
        if (lengthInInts > data.length) {
            int[] newData = new int[data.length * 2 + Math.max(0, lengthInInts - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }
}
