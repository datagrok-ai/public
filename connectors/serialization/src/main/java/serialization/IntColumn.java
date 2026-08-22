package serialization;

import serialization.codecs.BitIntList;
import serialization.codecs.IntRle;
import serialization.codecs.IntSequencePattern;

public class IntColumn extends AbstractColumn<Integer> {
    private static final String TYPE = Types.INT;
    public static final int None = -2147483648;

    // Best-encoder selection (ids 2/3/4) as in Dart; false restores int:raw only.
    public static boolean ADVANCED_ENCODERS =
            Boolean.parseBoolean(System.getProperty("grok.connect.advancedEncoders", "true"));

    private int[] data;

    public IntColumn(String name) {
        super(name);
        data = new int[initColumnSize];
    }

    public IntColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        data = new int[initColumnSize];
    }

    public IntColumn(String name, Integer[] values) {
        super(name);
        data = new int[initColumnSize];
        addAll(values);
    }

    public IntColumn(String name, int[] values) {
        super(name);
        data = values;
        length = values.length;
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
        Encoding.choose(data, length).write(buf, data, length);
    }

    // The encoder encode() picks for a value list and the payload bytes it takes after the 4-byte id
    // (codec headers included), so nested writers (DateTimeColumn components) can estimate before committing.
    public static final class Encoding {
        private static final int LIST_HEADER = 2 + 8;
        private static final int BIT_LIST_HEADER = 4 + 8 + 8 + 1 + LIST_HEADER;

        public final int id;
        public final int size;
        private final IntSequencePattern pattern;
        private final IntRle rle;
        private final int min;

        private Encoding(int id, int size, IntSequencePattern pattern, IntRle rle, int min) {
            this.id = id;
            this.size = size;
            this.pattern = pattern;
            this.rle = rle;
            this.min = min;
        }

        private static Encoding raw(int length) {
            return new Encoding(1, 1 + LIST_HEADER + length * 4, null, null, 0);
        }

        // Dart's getBestEncodingEstimate over IntMeta.encoders (pattern, raw, rle, bitIntList): skip -1, stop at 0, else smallest.
        public static Encoding choose(int[] data, int length) {
            if (!ADVANCED_ENCODERS)
                return raw(length);

            IntSequencePattern pattern = IntSequencePattern.fromList(data, length);
            if (pattern != null)
                return new Encoding(2, 5 * 4, pattern, null, 0);

            long min = Long.MAX_VALUE;
            long max = Long.MIN_VALUE;
            for (int i = 0; i < length; i++) {
                if (data[i] == None)
                    continue;
                min = Math.min(min, data[i]);
                max = Math.max(max, data[i]);
            }
            if (max < min)
                return raw(length);

            int best = length * 4;
            int bitListSize = BitIntList.sizeInBytes(length, min, max);
            // Range >= 2^30: bitIntList cannot beat raw and rle's deltas (up to twice the range) rarely fit its 31 bits, so skip its O(n) estimate.
            IntRle rle = BitIntList.msb(max - min) >= 31 ? null : new IntRle();
            int rleSize = rle == null ? -1 : rle.estimate(data, length, (int) min);
            if (rleSize != -1 && rleSize < best && rleSize <= bitListSize)
                return new Encoding(3, 3 * 4 + 2 * BIT_LIST_HEADER + rleSize, null, rle, (int) min);
            if (bitListSize < best)
                return new Encoding(4, BIT_LIST_HEADER + bitListSize, null, null, 0);
            return raw(length);
        }

        public void write(BufferAccessor buf, int[] data, int length) {
            buf.writeInt32(id);
            switch (id) {
                case 2:
                    pattern.serialize(buf);
                    break;
                case 3:
                    rle.encode(buf, data, length, min);
                    break;
                case 4:
                    BitIntList.fromList(data, 0, length).serialize(buf);
                    break;
                default:
                    buf.writeInt8((byte) 0);
                    buf.writeInt32List(data, 0, length);
            }
        }
    }

    @Override
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        switch (id) {
            case 1: // raw
                if (buf.readInt8() == ColumnEncoderArchiveType.ARCHIVE_TYPE_ZLIB)
                    data = ByteData.toInt32List(Zlib.inflate(buf.readUint8List()));
                else
                    data = buf.readInt32List();
                break;
            case 2: // pattern
                data = new serialization.codecs.IntSequencePattern(buf).toInt32List();
                break;
            case 3: // rle
                data = serialization.codecs.IntRle.decode(buf);
                break;
            case 4: // bitIntList
                data = serialization.codecs.BitIntList.fromBuffer(buf).toInt32List();
                break;
            default:
                throw new RuntimeException("decoding " + name + ": int encoder " + id + " not found");
        }
        length = data.length;
    }

    @Override
    public void add(Integer value) {
        ensureSpace(1);
        data[length++] = (value != null) ? value : None;
    }

    public void add(int value) {
        ensureSpace(1);
        data[length++] = value;
    }

    @Override
    public void addAll(Integer[] values) {
        ensureSpace(values.length);
        for (Integer value : values)
            data[length++] = (value != null) ? value : None;
    }

    @Override
    public Integer get(int idx) {
        return data[idx];
    }

    @Override
    public void set(int index, Integer value) {
        data[index] = value;
    }

    @Override
    public long memoryInBytes() {
        return (long) data.length * 4;
    }

    @Override
    public boolean isNone(int idx) {
        return data[idx] == None;
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            int[] newData = new int[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    @Override
    public Object toArray() {
        return data;
    }
}
