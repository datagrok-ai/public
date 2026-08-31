package serialization;

// Long storage until the first non-integral value (then decimal strings, today's wire format).
// With downcastAllowed the type is `int` while every value fits int32; sticky survives empty().
public class BigIntColumn extends AbstractColumn<Object> {
    static final String TYPE = Types.BIG_INT;
    public static final long None = Long.MIN_VALUE;

    private long[] longs;
    private long min = Long.MAX_VALUE;
    private long max = Long.MIN_VALUE;
    private StringColumn strings;
    private boolean downcastAllowed;
    private boolean sticky;

    public BigIntColumn(String name) {
        this(name, DEFAULT_INIT_SIZE);
    }

    public BigIntColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        longs = new long[initColumnSize];
    }

    public BigIntColumn(String name, String[] values) {
        this(name);
        addAll(values);
    }

    public void setDowncastAllowed(boolean value) {
        downcastAllowed = value;
    }

    public boolean isSticky() {
        return sticky;
    }

    public void setSticky(boolean value) {
        sticky = value;
    }

    @Override
    public String getType() {
        return downcastAllowed && !sticky && strings == null && length > 0 && max >= min
                && min > Integer.MIN_VALUE && max <= Integer.MAX_VALUE ? Types.INT : TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        longs = new long[initColumnSize];
        strings = null;
        min = Long.MAX_VALUE;
        max = Long.MIN_VALUE;
    }

    @Override
    public void encode(BufferAccessor buf) {
        if (getType().equals(Types.INT)) {
            IntColumn ints = new IntColumn(name, length);
            for (int i = 0; i < length; i++)
                ints.add(longs[i] == None ? IntColumn.None : (int) longs[i]);
            ints.encode(buf);
            return;
        }
        buf.writeInt32(1);
        toStrings().encode(buf);
    }

    @Override
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        StringColumn decoded = new StringColumn(name, 0);
        switch (id) {
            case 1: // bigInt:raw delegates to the nested string-column (categories) decode.
                decoded.decode(buf);
                break;
            case 2: // bigInt:list - polynomial {sizes, values, signs}.
                adoptFlat(decoded, decodeList(buf));
                break;
            case 3: // bigInt:capped - fixed-length polynomial columns.
                adoptFlat(decoded, decodeCapped(buf));
                break;
            default:
                throw new RuntimeException("decoding " + name + ": bigint encoder " + id + " not supported");
        }
        strings = decoded;
        longs = null;
        length = decoded.getLength();
    }

    @Override
    public void add(Object value) {
        if (strings == null) {
            if (value == null) {
                addLong(None);
                return;
            }
            if (isExactIntegral(value) && ((Number) value).longValue() != None) {
                addLong(((Number) value).longValue());
                return;
            }
            switchToStrings();
        }
        strings.add(toText(value));
        length++;
    }

    @Override
    public void addAll(Object[] values) {
        for (Object value : values)
            add(value);
    }

    @Override
    public Object get(int idx) {
        if (strings != null)
            return strings.get(idx);
        return longs[idx] == None ? null : (Object) longs[idx];
    }

    @Override
    public void set(int index, Object value) {
        if (strings == null && (value == null || isExactIntegral(value) && ((Number) value).longValue() != None)) {
            longs[index] = value == null ? None : ((Number) value).longValue();
            track(longs[index]);
            return;
        }
        if (strings == null)
            switchToStrings();
        strings.set(index, toText(value));
    }

    @Override
    public long memoryInBytes() {
        return strings != null ? strings.memoryInBytes() : (long) longs.length * 8;
    }

    @Override
    public boolean isNone(int idx) {
        return strings != null ? strings.isNone(idx) : longs[idx] == None;
    }

    @Override
    public Object toArray() {
        return strings != null ? strings.toArray() : longs;
    }

    private static boolean isExactIntegral(Object value) {
        return value instanceof Long || value instanceof Integer || value instanceof Short || value instanceof Byte;
    }

    private static String toText(Object value) {
        return value == null ? null : value.toString();
    }

    private void addLong(long value) {
        if (length >= longs.length) {
            long[] newData = new long[Math.max(longs.length * 2, length + 1)];
            System.arraycopy(longs, 0, newData, 0, longs.length);
            longs = newData;
        }
        longs[length++] = value;
        track(value);
    }

    private void track(long value) {
        if (value == None)
            return;
        if (value < min) min = value;
        if (value > max) max = value;
        if (value <= Integer.MIN_VALUE || value > Integer.MAX_VALUE)
            sticky = true;
    }

    private StringColumn toStrings() {
        if (strings != null)
            return strings;
        StringColumn result = new StringColumn(name, length);
        for (int i = 0; i < length; i++)
            result.add(longs[i] == None ? null : Long.toString(longs[i]));
        return result;
    }

    private void switchToStrings() {
        strings = toStrings();
        longs = null;
        sticky = true;
    }

    // Stores flat decimal strings (null = None) with an identity index map, so get(idx)
    // returns values[idx] without re-categorizing.
    private static void adoptFlat(StringColumn target, String[] values) {
        int[] indices = new int[values.length];
        for (int i = 0; i < values.length; i++)
            indices[i] = i;
        target.adoptDecoded(values, indices);
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
