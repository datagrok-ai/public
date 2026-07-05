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

    @Override
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        switch (id) {
            case 1: // bigInt:raw delegates to the nested string-column (categories) decode.
                super.decode(buf);
                break;
            case 2: // bigInt:list - polynomial {sizes, values, signs}.
                adoptFlat(decodeList(buf));
                break;
            case 3: // bigInt:capped - fixed-length polynomial columns.
                adoptFlat(decodeCapped(buf));
                break;
            default:
                throw new RuntimeException("decoding " + name + ": bigint encoder " + id + " not supported");
        }
    }

    // Stores flat decimal strings (null = None) into the StringColumn parent with an
    // identity index map, so get(idx) returns values[idx] without re-categorizing.
    private void adoptFlat(String[] values) {
        int[] indices = new int[values.length];
        for (int i = 0; i < values.length; i++)
            indices[i] = i;
        adoptDecoded(values, indices);
    }

    // Ports BigIntListEncoder.decode (big_int_column_encoders.dart:51-60).
    private static String[] decodeList(BufferAccessor buf) {
        Column<?> sizesCol = buf.readColumn();
        Column<?> valuesCol = buf.readColumn();
        Column<?> signsCol = buf.readColumn();
        int n = sizesCol.getLength();
        int[] sizes = (int[]) sizesCol.toArray();
        int[] valData = (int[]) valuesCol.toArray();
        String[] out = new String[n];
        int vi = 0;
        for (int i = 0; i < n; i++) {
            int sz = sizes[i];
            if (sz > 0) {
                int[] coeffs = new int[sz];
                for (int k = 0; k < sz; k++)
                    coeffs[k] = valData[vi++];
                out[i] = bigIntToString(coeffs, (Boolean) signsCol.get(i));
            }
        }
        return out;
    }

    // Ports BigIntCappedEncoder.decode (big_int_column_encoders.dart:111-130).
    private static String[] decodeCapped(BufferAccessor buf) {
        int bigIntLength = buf.readUint16();
        int columnLength = buf.readInt32();
        int batchSize = 30000;
        String[] out = new String[columnLength];
        for (int start = 0; start < columnLength; start += batchSize) {
            if (start + 2 * batchSize >= columnLength)
                batchSize = columnLength - start;
            Column<?>[] valuesColumns = new Column<?>[bigIntLength];
            for (int i = 0; i < bigIntLength; i++)
                valuesColumns[i] = buf.readColumn();
            for (int nn = 0; nn < batchSize; nn++) {
                int[] values = new int[bigIntLength];
                for (int i = 0; i < bigIntLength; i++)
                    values[i] = (Integer) valuesColumns[i].get(nn);
                if (values[0] != -1)
                    out[start + nn] = bigIntToString(values, true);
            }
        }
        return out;
    }

    // Reconstructs the decimal string from little-endian base-10^7 coefficients
    // (ddt/lib/src/utils/big_int.dart BigInt.toString, _LOG_BASE=7, _BASE=10^7).
    // Trailing (most-significant) zero coefficients are dropped as the BigInt
    // constructor does.
    static String bigIntToString(int[] coeffs, boolean isPositive) {
        int len = coeffs.length;
        while (len > 0 && coeffs[len - 1] == 0)
            len--;
        if (len == 0)
            return "0";
        StringBuilder sb = new StringBuilder();
        if (!isPositive)
            sb.append('-');
        sb.append(coeffs[len - 1]);
        for (int i = len - 2; i >= 0; i--)
            sb.append(String.format("%07d", coeffs[i]));
        return sb.toString();
    }
}
