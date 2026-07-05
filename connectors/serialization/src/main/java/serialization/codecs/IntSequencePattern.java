package serialization.codecs;

import serialization.BufferAccessor;

// Decode-only port of ddt/lib/src/serialization/int_pattern.dart:40-59.
// Reconstructs data that consists of repeating blocks with varying periodicity.
public class IntSequencePattern {
    private int start;
    private int step;
    private int blockLength;
    private int blocksPerCycle;
    private int length;

    public IntSequencePattern(BufferAccessor buffer) {
        start = buffer.readInt32();
        step = buffer.readInt32();
        blockLength = buffer.readInt32();
        blocksPerCycle = buffer.readInt32();
        length = buffer.readInt32();
    }

    public int get(int i) {
        return start + ((i % (blocksPerCycle * blockLength)) / blockLength) * step;
    }

    public int[] toInt32List() {
        int[] r = new int[length];
        for (int i = 0; i < length; i++)
            r[i] = get(i);
        return r;
    }
}
