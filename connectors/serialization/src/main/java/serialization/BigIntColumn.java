package serialization;

// Big integer column.
public class BigIntColumn extends StringColumn {
    static final String TYPE = Types.BIG_INT;

    public BigIntColumn() {
        super();
    }

    public BigIntColumn(int initColumnSize) {
        super(initColumnSize);
    }

    public BigIntColumn(String[] values) {
        super(values);
    }

    public String getType() {
        return TYPE;
    }

    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);
        super.encode(buf);
    }
}
