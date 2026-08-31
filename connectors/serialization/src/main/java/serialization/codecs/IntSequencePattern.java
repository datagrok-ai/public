package serialization.codecs;

import serialization.BufferAccessor;

// Port of ddt/lib/src/serialization/int_pattern.dart.
// Data that consists of repeating blocks with varying periodicity, e.g.
// [4, 5, 6, 7, 8], [5, 7, 9, 11], [1, 1, 2, 2, 3, 3, 4, 4].
public class IntSequencePattern {
    private int start;
    private long step;
    private int blockLength;
    private int blocksPerCycle;
    private int length;

    private IntSequencePattern(int start, long step, int blockLength, int blocksPerCycle, int length) {
        this.start = start;
        this.step = step;
        this.blockLength = blockLength;
        this.blocksPerCycle = blocksPerCycle;
        this.length = length;
    }

    public IntSequencePattern(BufferAccessor buffer) {
        start = buffer.readInt32();
        step = buffer.readInt32();
        blockLength = buffer.readInt32();
        blocksPerCycle = buffer.readInt32();
        length = buffer.readInt32();
    }

    // Deltas are kept 64-bit like Dart ints, so a wrapped int32 difference never matches.
    public static IntSequencePattern fromList(int[] values, int length) {
        if (length < 3)
            return null;

        int start = values[0];
        long step = (long) values[1] - values[0];
        int blocksPerCycle = -1;
        int blockLength = (step == 0 ? -1 : 1);

        for (int i = 2; i < length; i++) {
            int x = values[i];

            if (blockLength == -1) {
                if (x != start) {
                    step = (long) x - start;
                    blockLength = i;
                }
                continue;
            }

            if (blocksPerCycle > 0) {
                long expected = start + ((i % (blocksPerCycle * blockLength)) / blockLength) * step;
                if (x == expected)
                    continue;
                else
                    return null;
            }

            boolean blockStart = (i % blockLength) == 0;
            if (!blockStart) {
                if (x == values[i - 1])
                    continue;
                else
                    return null;
            }
            else {
                if (x == start) {
                    blocksPerCycle = i / blockLength;
                    continue;
                }
                if (x == values[i - 1] + step)
                    continue;
            }

            return null;
        }

        return new IntSequencePattern(start, step, blockLength, blocksPerCycle == -1 ? length : blocksPerCycle, length);
    }

    public void serialize(BufferAccessor buffer) {
        buffer.writeInt32(start);
        buffer.writeInt32((int) step);
        buffer.writeInt32(blockLength);
        buffer.writeInt32(blocksPerCycle);
        buffer.writeInt32(length);
    }

    public int get(int i) {
        return (int) (start + ((i % ((long) blocksPerCycle * blockLength)) / blockLength) * step);
    }

    public int[] toInt32List() {
        int[] r = new int[length];
        for (int i = 0; i < length; i++)
            r[i] = get(i);
        return r;
    }
}
