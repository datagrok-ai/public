package serialization;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import java.util.Arrays;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

// Benchmark results (avg ms, 10 measured runs after 5 warmup runs):
//
// Scenario                        Old (sort-based)    New (HashMap-based)    Speedup
// 100K rows / 10 categories            10.2ms                1.2ms            8.5x
// 100K rows / 100 categories           13.4ms                1.6ms            8.4x
// 100K rows / 1K categories            18.4ms                1.6ms           11.5x
// 500K rows / 10 categories            46.5ms                3.8ms           12.2x
// 500K rows / 100 categories           74.6ms                4.3ms           17.3x
// 500K rows / 1K categories           100.7ms                7.4ms           13.6x
// 1M rows / 10 categories              98.8ms                7.6ms           13.0x
// 1M rows / 100 categories            154.6ms                8.3ms           18.6x
// 1M rows / 1K categories             212.7ms               14.6ms           14.6x
// 1M rows / 10K categories            304.1ms               18.0ms           16.9x
// 100K rows / all unique               31.1ms               16.3ms            1.9x
// 500K rows / all unique              269.4ms              110.2ms            2.4x
// 1M rows / all unique                814.4ms              247.8ms            3.3x
public class StringColumnPerformanceTest {
    private static final int WARMUP_RUNS = 5;
    private static final int MEASURED_RUNS = 10;

    @ParameterizedTest(name = "rows={0}, categories={1}")
    @CsvSource({
        "100000,   10,   '100K rows / 10 categories'",
        "100000,   100,  '100K rows / 100 categories'",
        "100000,   1000, '100K rows / 1K categories'",
        "500000,   10,   '500K rows / 10 categories'",
        "500000,   100,  '500K rows / 100 categories'",
        "500000,   1000, '500K rows / 1K categories'",
        "1000000,  10,   '1M rows / 10 categories'",
        "1000000,  100,  '1M rows / 100 categories'",
        "1000000,  1000, '1M rows / 1K categories'",
        "1000000,  10000, '1M rows / 10K categories'",
        "100000,   100000,'100K rows / all unique'",
        "500000,   500000,'500K rows / all unique'",
        "1000000,  1000000,'1M rows / all unique'",
    })
    void encodeBenchmark(int rowCount, int categoryCount, String label) {
        StringColumn col = buildColumn(rowCount, categoryCount);

        // Warmup
        for (int i = 0; i < WARMUP_RUNS; i++)
            encodeColumn(col);

        // Measured runs
        long[] times = new long[MEASURED_RUNS];
        for (int i = 0; i < MEASURED_RUNS; i++) {
            long start = System.nanoTime();
            encodeColumn(col);
            times[i] = System.nanoTime() - start;
        }

        long min = Long.MAX_VALUE, max = 0, sum = 0;
        for (long t : times) {
            if (t < min) min = t;
            if (t > max) max = t;
            sum += t;
        }
        double avgMs = (sum / (double) MEASURED_RUNS) / 1_000_000.0;
        double minMs = min / 1_000_000.0;
        double maxMs = max / 1_000_000.0;

        System.out.printf("%-30s  avg=%.1fms  min=%.1fms  max=%.1fms%n",
                label, avgMs, minMs, maxMs);
    }

    @Test
    void encodeCorrectnessAfterCategorize() {
        StringColumn col = new StringColumn("test", new String[]{"banana", "apple", null, "banana", "", "apple", "cherry"});
        byte[] encoded = encodeColumn(col);
        assertTrue(encoded.length > 0, "Encoded output should not be empty");

        // Encode again to verify deterministic output
        byte[] encoded2 = encodeColumn(col);
        assertArrayEquals(encoded, encoded2, "Encoding should be deterministic");
    }

    @Test
    void encodeSingleCategory() {
        String[] values = new String[10000];
        Arrays.fill(values, "same");
        StringColumn col = new StringColumn("single", values);
        byte[] encoded = encodeColumn(col);
        assertTrue(encoded.length > 0);
    }

    @Test
    void encodeAllUnique() {
        String[] values = new String[10000];
        for (int i = 0; i < values.length; i++)
            values[i] = "val_" + i;
        StringColumn col = new StringColumn("unique", values);
        byte[] encoded = encodeColumn(col);
        assertTrue(encoded.length > 0);
    }

    @Test
    void encodeWithNulls() {
        String[] values = new String[10000];
        Random rng = new Random(42);
        for (int i = 0; i < values.length; i++)
            values[i] = rng.nextBoolean() ? null : "cat_" + (i % 5);
        StringColumn col = new StringColumn("nulls", values);
        byte[] encoded = encodeColumn(col);
        assertTrue(encoded.length > 0);
    }

    private static StringColumn buildColumn(int rowCount, int categoryCount) {
        String[] categories = new String[categoryCount];
        for (int i = 0; i < categoryCount; i++)
            categories[i] = "category_" + i;

        StringColumn col = new StringColumn("bench", rowCount);
        Random rng = new Random(12345);
        for (int i = 0; i < rowCount; i++)
            col.add(categories[rng.nextInt(categoryCount)]);
        return col;
    }

    private static byte[] encodeColumn(StringColumn col) {
        BufferAccessor buf = new BufferAccessor();
        col.encode(buf);
        return buf.toUint8List();
    }
}
