package serialization;

import serialization.codecs.BitIntList;
import serialization.codecs.BitPacking;
import serialization.codecs.IntRle;
import serialization.codecs.IntSequencePattern;

public class IntColumn extends AbstractColumn<Integer> {
    private static final String TYPE = Types.INT;
    public static final int None = -2147483648;

    // Best-encoder selection (ids 2/3/4, level-2 ids 5/6 - GROK-20761) as in Dart; false restores int:raw only.
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
        private long[] packMinMax;   // id 5: [min, max, hasNone]
        private long[] deltaRange;   // id 6: [dMin, dMax, 0]
        private int none;            // id 6: the None substitute (min - 1)

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
            long hasNone = 0;
            for (int i = 0; i < length; i++) {
                if (data[i] == None) {
                    hasNone = 1;
                    continue;
                }
                min = Math.min(min, data[i]);
                max = Math.max(max, data[i]);
            }
            long[] packMinMax = max < min ? null : new long[]{min, max, hasNone};

            // Dart iteration order (raw, rle, bitIntList, bitPacked, delta), strict-min wins.
            int bestSize = length * 4;
            int bestId = 1;
            int bitListSize = max < min ? Integer.MAX_VALUE : BitIntList.sizeInBytes(length, min, max);
            // Range >= 2^30: bitIntList cannot beat raw and rle's deltas (up to twice the range) rarely fit its 31 bits, so skip its O(n) estimate.
            IntRle rle = max < min || BitIntList.msb(max - min) >= 31 ? null : new IntRle();
            int rleSize = rle == null ? -1 : rle.estimate(data, length, (int) min);
            if (rleSize != -1 && rleSize < bestSize) {
                bestId = 3;
                bestSize = rleSize;
            }
            if (bitListSize < bestSize) {
                bestId = 4;
                bestSize = bitListSize;
            }
            int bpBits = BitPacking.bitsFor(packMinMax);
            if (bpBits >= 0 && BitPacking.sizeInBytes(length, bpBits) < bestSize) {
                bestId = 5;
                bestSize = BitPacking.sizeInBytes(length, bpBits);
            }
            long[] deltaRange = null;
            int none = packMinMax == null ? 0 : (int) (min - 1);
            if (length >= 2) {
                long dMin = 0, dMax = 0;
                boolean any = false;
                boolean fits = true;
                long prev = data[0] == None ? none : data[0];
                for (int i = 1; i < length && fits; i++) {
                    long v = data[i] == None ? none : data[i];
                    long d = v - prev;
                    fits = d >= Integer.MIN_VALUE && d <= Integer.MAX_VALUE;
                    prev = v;
                    dMin = any ? Math.min(dMin, d) : d;
                    dMax = any ? Math.max(dMax, d) : d;
                    any = true;
                }
                long[] range = new long[]{dMin, dMax, 0};
                int dBits = fits ? BitPacking.bitsFor(range) : -1;
                if (dBits >= 0 && 8 + BitPacking.sizeInBytes(length - 1, dBits) < bestSize) {
                    bestId = 6;
                    bestSize = 8 + BitPacking.sizeInBytes(length - 1, dBits);
                    deltaRange = range;
                }
            }

            switch (bestId) {
                case 3:
                    return new Encoding(3, 3 * 4 + 2 * BIT_LIST_HEADER + rleSize, null, rle, (int) min);
                case 4:
                    return new Encoding(4, BIT_LIST_HEADER + bitListSize, null, null, 0);
                case 5: {
                    Encoding e = new Encoding(5, bestSize, null, null, 0);
                    e.packMinMax = packMinMax;
                    return e;
                }
                case 6: {
                    Encoding e = new Encoding(6, bestSize, null, null, 0);
                    e.deltaRange = deltaRange;
                    e.none = none;
                    return e;
                }
                default:
                    return raw(length);
            }
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
                case 5:
                    BitPacking.write(buf, data, 0, length, packMinMax);
                    break;
                case 6: {
                    int first = data[0] == None ? none : data[0];
                    buf.writeInt32(first);
                    buf.writeInt32(none);
                    int[] deltas = new int[length - 1];
                    int prev = first;
                    for (int i = 1; i < length; i++) {
                        int v = data[i] == None ? none : data[i];
                        deltas[i - 1] = v - prev;
                        prev = v;
                    }
                    BitPacking.write(buf, deltas, 0, deltas.length, deltaRange);
                    break;
                }
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
            case 5: // bitPacked (dense FOR packing, GROK-20761)
                data = BitPacking.read(buf);
                break;
            case 6: { // delta (first + bit-packed deltas, GROK-20761)
                int first = buf.readInt32();
                int none = buf.readInt32();
                int[] deltas = BitPacking.read(buf);
                data = new int[deltas.length + 1];
                int acc = first;
                data[0] = acc == none ? None : acc;
                for (int i = 0; i < deltas.length; i++) {
                    acc += deltas[i];
                    data[i + 1] = acc == none ? None : acc;
                }
                break;
            }
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
