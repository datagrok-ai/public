package serialization;

public class BoolColumn extends AbstractColumn<Boolean> {
    private static final String TYPE = Types.BOOL;

    private int[] data;

    public BoolColumn(String name) {
        super(name);
        data = new int[initColumnSize];
    }

    public BoolColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        data = new int[initColumnSize];
    }

    public BoolColumn(String name, Boolean[] values) {
        super(name);
        data = new int[initColumnSize];
        addAll(values);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        data = new int[initColumnSize];
    }

    @Override
    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);
        buf.writeInt64(length);
        buf.writeInt8((byte)0);
        buf.writeUint32List(data, 0, ((length + 0x1F) / 0x20));
    }

    @Override
    public void add(Boolean value) {
        ensureSpace(1);
        if ((value != null) && value)
            data[length / 0x20] |= 1 << ((length % 0x20) & 0x1F);
        length++;
    }

    @Override
    public void addAll(Boolean[] values) {
        ensureSpace(values.length);
        for (Boolean value : values) {
            if ((value != null) && value)
                data[length / 0x20] |= 1 << ((length % 0x20) & 0x1F);
            length++;
        }
    }

    @Override
    public Boolean get(int idx) {
        return (data[idx / 0x20] & (1 << (idx % 0x20 & 0x1F))) != 0;
    }

    @Override
    public void set(int index, Boolean value) {
        if (value != null && value)
            data[index / 0x20] |= 1 << (index % 0x20 & 0x1F);
        else
            data[index / 0x20] &= ~(1 << (index % 0x20 & 0x1F));
    }

    @Override
    public boolean isNone(int idx) {
        return false;
    }

    @Override
    public long memoryInBytes() {
        return (long) data.length * 4;
    }

    private void ensureSpace(int extraLength) {
        int lengthInInts = ((length + extraLength + 0x1F) / 0x20);
        if (lengthInInts > data.length) {
            int[] newData = new int[data.length * 2 + Math.max(0, lengthInInts - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    @Override
    public Object toArray() {
        return data;
    }
}
