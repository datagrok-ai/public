package serialization;

import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;


public class BufferAccessorPerformanceTest {
    private static final int WARMUP_RUNS = 5;
    private static final int MEASURED_RUNS = 10;

    @Test
    void writeRawBytesBenchmark() {
        byte[] data = new byte[10_000_000];
        new Random(42).nextBytes(data);

        for (int i = 0; i < WARMUP_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor();
            buf.writeRawBytes(data);
        }

        long[] times = new long[MEASURED_RUNS];
        for (int i = 0; i < MEASURED_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor();
            long start = System.nanoTime();
            buf.writeRawBytes(data);
            times[i] = System.nanoTime() - start;
        }
        printTimes("writeRawBytes 10MB", times);
    }

    @Test
    void writeUint8ListBenchmark() {
        byte[] data = new byte[10_000_000];
        new Random(42).nextBytes(data);

        for (int i = 0; i < WARMUP_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor();
            buf.writeUint8List(data);
        }

        long[] times = new long[MEASURED_RUNS];
        for (int i = 0; i < MEASURED_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor();
            long start = System.nanoTime();
            buf.writeUint8List(data);
            times[i] = System.nanoTime() - start;
        }
        printTimes("writeUint8List 10MB", times);
    }

    @Test
    void writeStringListBenchmark() {
        String[] strings = new String[100_000];
        Random rng = new Random(42);
        for (int i = 0; i < strings.length; i++) {
            StringBuilder sb = new StringBuilder(50);
            for (int j = 0; j < 50; j++)
                sb.append((char) ('a' + rng.nextInt(26)));
            strings[i] = sb.toString();
        }

        for (int i = 0; i < WARMUP_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor();
            buf.writeStringList(strings);
        }

        long[] times = new long[MEASURED_RUNS];
        for (int i = 0; i < MEASURED_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor();
            long start = System.nanoTime();
            buf.writeStringList(strings);
            times[i] = System.nanoTime() - start;
        }
        printTimes("writeStringList 100K×50B", times);
    }

    @Test
    void insertStringHeaderBenchmark() {
        // Simulate a typical streamed DataFrame chunk (~5MB)
        byte[] blob = new byte[5_000_000];
        new Random(42).nextBytes(blob);
        String header = "{\"type\":\"DataQueryRunResult\",\"log\":\"\",\"columns\":10,\"rows\":1000000}";

        for (int i = 0; i < WARMUP_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor(blob.clone());
            buf.bufPos = blob.length;
            buf.insertStringHeader(header);
            buf.toUint8List();
        }

        long[] times = new long[MEASURED_RUNS];
        for (int i = 0; i < MEASURED_RUNS; i++) {
            BufferAccessor buf = new BufferAccessor(blob.clone());
            buf.bufPos = blob.length;
            long start = System.nanoTime();
            buf.insertStringHeader(header);
            buf.toUint8List();
            times[i] = System.nanoTime() - start;
        }
        printTimes("insertHeader+toBytes 5MB", times);
    }

    @Test
    void fullDataFrameBenchmark() {
        DataFrame df = new DataFrame();
        int rowCount = 500_000;

        IntColumn intCol = new IntColumn("ints", rowCount);
        FloatColumn floatCol = new FloatColumn("floats", rowCount);
        StringColumn strCol = new StringColumn("strings", rowCount);
        Random rng = new Random(42);
        String[] cats = new String[100];
        for (int i = 0; i < cats.length; i++)
            cats[i] = "category_" + i;
        for (int i = 0; i < rowCount; i++) {
            intCol.add(rng.nextInt());
            floatCol.add(rng.nextFloat());
            strCol.add(cats[rng.nextInt(cats.length)]);
        }
        df.addColumns(new Column[]{intCol, floatCol, strCol});

        for (int i = 0; i < WARMUP_RUNS; i++)
            df.toByteArray();

        long[] times = new long[MEASURED_RUNS];
        for (int i = 0; i < MEASURED_RUNS; i++) {
            long start = System.nanoTime();
            df.toByteArray();
            times[i] = System.nanoTime() - start;
        }
        printTimes("DataFrame 500K×3cols", times);
    }

    @Test
    void correctnessWriteRawBytes() {
        byte[] data = {1, 2, 3, 4, 5, 127, -128, 0};
        BufferAccessor buf = new BufferAccessor();
        buf.writeRawBytes(data);
        byte[] result = buf.toUint8List();
        assertArrayEquals(data, result);
    }

    @Test
    void correctnessInsertStringHeader() {
        BufferAccessor buf = new BufferAccessor();
        buf.writeInt32(42);
        byte[] before = buf.toUint8List();

        buf.insertStringHeader("hello");
        byte[] after = buf.toUint8List();

        // The last 4 bytes should still be the int32(42)
        byte[] tail = new byte[4];
        System.arraycopy(after, after.length - 4, tail, 0, 4);
        assertArrayEquals(before, tail, "Original data must be preserved after header insert");
        assertTrue(after.length > before.length, "Buffer should grow after header insert");
    }

    private static void printTimes(String label, long[] times) {
        long min = Long.MAX_VALUE, max = 0, sum = 0;
        for (long t : times) {
            if (t < min) min = t;
            if (t > max) max = t;
            sum += t;
        }
        double avgMs = (sum / (double) times.length) / 1_000_000.0;
        double minMs = min / 1_000_000.0;
        double maxMs = max / 1_000_000.0;
        System.out.printf("%-35s  avg=%.1fms  min=%.1fms  max=%.1fms%n",
                label, avgMs, minMs, maxMs);
    }
}
