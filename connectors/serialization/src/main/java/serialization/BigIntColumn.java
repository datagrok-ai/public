package serialization;


// Big integer column.
public class BigIntColumn extends StringColumn {
    static final String TYPE = Types.BIG_INT;

    public String getType() {
        return TYPE;
    }

    public BigIntColumn() {
        super();
    }

    public BigIntColumn(String[] values) {
        super(values);
    }

    public void encode(BufferAccessor buf) {
        buf.writeInt32(1);
        super.encode(buf);
    }
}
