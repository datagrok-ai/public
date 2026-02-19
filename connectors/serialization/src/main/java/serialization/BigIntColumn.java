package serialization;

public class BigIntColumn extends StringColumn {
    static final String TYPE = Types.BIG_INT;

    public BigIntColumn(String name) {
        super(name);
    }

    public BigIntColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
    }

    public BigIntColumn(String name, String[] values) {
        super(name, values);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);
        super.encode(buf);
    }
}
