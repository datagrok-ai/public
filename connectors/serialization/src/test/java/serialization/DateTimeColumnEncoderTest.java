package serialization;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.util.Random;

// DateTimeColumn's dateTime:int (id 2) writer: byte-exact round trip through the Java id 2 reader
// (decodeComponents) and the selection rule against dateTime:microseconds (id 3).
public class DateTimeColumnEncoderTest {
    private static final long START = micros(2023, 1, 1, 8, 0, 0, 0);
    private static final long YEAR_0001 = micros(1, 1, 1, 0, 0, 0, 0);

    @AfterEach
    public void restoreSwitches() {
        DateTimeColumn.COMPONENT_ENCODER = true;
        IntColumn.ADVANCED_ENCODERS = true;
    }

    static long micros(int y, int mo, int d, int h, int mi, int s, int micro) {
        return LocalDateTime.of(y, mo, d, h, mi, s, micro * 1000).toEpochSecond(ZoneOffset.UTC) * 1_000_000L + micro;
    }

    static Double[] regular(int rows) {
        Double[] v = new Double[rows];
        for (int i = 0; i < rows; i++)
            v[i] = (double) (START + i * 1_500_000L);
        return v;
    }

    static Double[] jitter(int rows) {
        Random rnd = new Random(42);
        Double[] v = regular(rows);
        for (int i = 0; i < rows; i++)
            v[i] += rnd.nextInt(5001) * 1000L;
        return v;
    }

    static Double[] dateOnly(int rows) {
        Double[] v = new Double[rows];
        for (int i = 0; i < rows; i++)
            v[i] = (double) micros(2023, 1, 1, 0, 0, 0, 0) + (i / 4) * 86_400_000_000L;
        return v;
    }

    static Double[] mixedNone(int rows) {
        Double[] v = regular(rows);
        for (int i = 0; i < rows; i += 7)
            v[i] = null;
        return v;
    }

    static Double[] subSecond(int rows) {
        Double[] v = regular(rows);
        for (int i = 0; i < rows; i++)
            v[i] += (i * 7) % 1000;
        return v;
    }

    // Whole-millisecond timestamps spread uniformly over 1970-2100.
    static Double[] randomMillis(int rows) {
        Random rnd = new Random(7);
        Double[] v = new Double[rows];
        for (int i = 0; i < rows; i++)
            v[i] = (double) ((long) (rnd.nextDouble() * 130 * 365.25 * 86_400_000L) * 1000L);
        return v;
    }

    static Double[] randomMicros(int rows) {
        Random rnd = new Random(11);
        Double[] v = randomMillis(rows);
        for (int i = 0; i < rows; i++)
            v[i] += rnd.nextInt(1000);
        return v;
    }

    private static byte[] encode(DateTimeColumn col) {
        BufferAccessor buf = new BufferAccessor();
        col.encode(buf);
        return buf.toUint8List();
    }

    private static byte[] encodeComponents(DateTimeColumn col) {
        BufferAccessor buf = new BufferAccessor();
        Assertions.assertTrue(col.encodeComponents(buf, true));
        return buf.toUint8List();
    }

    private static DateTimeColumn decode(byte[] bytes) {
        DateTimeColumn col = new DateTimeColumn("", 0);
        col.decode(new BufferAccessor(bytes));
        return col;
    }

    private static void assertSame(Double[] expected, DateTimeColumn actual, String what) {
        Assertions.assertEquals(expected.length, actual.getLength(), what + " length");
        for (int r = 0; r < expected.length; r++) {
            Assertions.assertEquals(expected[r] == null, actual.isNone(r), what + "[" + r + "] none");
            if (expected[r] != null)
                Assertions.assertEquals(expected[r], actual.get(r), 0.0, what + "[" + r + "]");
        }
    }

    private static void assertComponentsRoundTrip(String what, Double[] values) {
        byte[] bytes = encodeComponents(new DateTimeColumn(what, values));
        Assertions.assertEquals(2, new BufferAccessor(bytes).readInt32(), what + " encoder id");
        assertSame(values, decode(bytes), what);
    }

    private static int assertSelects(String what, Double[] values, int expectedId) {
        byte[] bytes = encode(new DateTimeColumn(what, values));
        int id = new BufferAccessor(bytes).readInt32();
        String sizes = what + ": id " + id + ", payload " + (bytes.length - 4) + " B vs length*8 = " + values.length * 8 + " B";
        System.out.println(sizes);
        Assertions.assertEquals(expectedId, id, sizes);
        if (id == 2)
            Assertions.assertTrue(bytes.length - 4 < values.length * 8, sizes);
        assertSame(values, decode(bytes), what);
        return bytes.length - 4;
    }

    @Test
    public void componentsRoundTripThroughTheJavaReader() {
        assertComponentsRoundTrip("regular", regular(5000));
        assertComponentsRoundTrip("jitter", jitter(5000));
        assertComponentsRoundTrip("dateOnly", dateOnly(5000));
        assertComponentsRoundTrip("allNone", new Double[1000]);
        assertComponentsRoundTrip("mixedNone", mixedNone(5000));
        assertComponentsRoundTrip("subSecond", subSecond(5000));
        assertComponentsRoundTrip("randomMicros", randomMicros(5000));
        assertComponentsRoundTrip("edges", new Double[] {
                -1.0, -1_000_000.0, 0.0, 1.0,
                (double) micros(1969, 12, 31, 23, 59, 59, 999_999),
                (double) micros(1900, 6, 15, 12, 34, 56, 789_000),
                (double) micros(9999, 12, 31, 23, 59, 59, 999_999),
                (double) micros(1, 1, 2, 0, 0, 0, 0),
                (double) micros(1, 1, 1, 0, 0, 0, 16),
                null});
    }

    @Test
    public void noneTupleIsIndistinguishableFromYearOne() {
        Double[] values = {(double) YEAR_0001, (double) START};
        DateTimeColumn decoded = decode(encodeComponents(new DateTimeColumn("y1", values)));
        Assertions.assertTrue(decoded.isNone(0), "0001-01-01T00:00:00Z is the None tuple (1,1,1,0,0,0,0,0)");
        Assertions.assertEquals(values[1], decoded.get(1), 0.0);
    }

    @Test
    public void selectsComponentsOnlyWhenSmallerAndLossless() {
        int regular = assertSelects("regular 1500 ms cadence", regular(100_000), 2);
        int jitter = assertSelects("jitter 0-5000 ms", jitter(100_000), 2);
        int dateOnly = assertSelects("date-only", dateOnly(100_000), 2);
        int randomMs = assertSelects("random ms 1970-2100", randomMillis(100_000), 2);
        Assertions.assertTrue(regular < jitter && jitter < randomMs, "regular " + regular + " < jitter " + jitter + " < random " + randomMs);
        Assertions.assertTrue(dateOnly < regular, "date-only " + dateOnly + " < regular " + regular);

        assertSelects("random us 1970-2100 (microseconds are lossy in dart2js)", randomMicros(100_000), 3);
        assertSelects("sub-second us", subSecond(100_000), 3);
        int mixedNone = assertSelects("mixed None (None rows break the cadence of every component)", mixedNone(100_000), 2);
        Assertions.assertTrue(randomMs > mixedNone && mixedNone > jitter, "random " + randomMs + " > mixed None " + mixedNone + " > jitter " + jitter);
        assertSelects("all None (constant (1,1,1) tuple)", new Double[100_000], 2);
        assertSelects("tiny column (headers dominate)", regular(6), 3);
        assertSelects("empty column", new Double[0], 3);
    }

    @Test
    public void switchesRestoreMicroseconds() {
        DateTimeColumn.COMPONENT_ENCODER = false;
        assertSelects("COMPONENT_ENCODER off", regular(10_000), 3);
        DateTimeColumn.COMPONENT_ENCODER = true;
        IntColumn.ADVANCED_ENCODERS = false;
        assertSelects("ADVANCED_ENCODERS off (raw components are 12+ B/row)", regular(10_000), 3);
    }

    @Test
    public void componentsMatchDartUtcBreakdown() {
        Double[] values = {(double) micros(2024, 2, 29, 23, 45, 6, 789_012), (double) micros(1969, 12, 31, 23, 59, 59, 999_999)};
        byte[] bytes = encodeComponents(new DateTimeColumn("c", values));
        BufferAccessor buf = new BufferAccessor(bytes);
        Assertions.assertEquals(2, buf.readInt32());
        int[][] expected = {{2024, 1969}, {2, 12}, {29, 31}, {23, 23}, {45, 59}, {6, 59}, {789, 999}, {12, 999}};
        for (int c = 0; c < 8; c++) {
            Assertions.assertEquals(1, buf.readInt8(), "component " + c + " present");
            IntColumn ic = new IntColumn("", 0);
            ic.decode(buf);
            Assertions.assertArrayEquals(expected[c], (int[]) ic.toArray(), "component " + c);
        }
        Assertions.assertEquals(bytes.length, buf.bufPos, "no trailing bytes");
    }

    @Test
    public void absentComponentsAreFlaggedZero() {
        byte[] bytes = encodeComponents(new DateTimeColumn("d", dateOnly(100)));
        BufferAccessor buf = new BufferAccessor(bytes);
        Assertions.assertEquals(2, buf.readInt32());
        for (int c = 0; c < 8; c++) {
            Assertions.assertEquals(c < 3 ? 1 : 0, buf.readInt8(), "component " + c + " present flag");
            if (c < 3)
                new IntColumn("", 0).decode(buf);
        }
        Assertions.assertEquals(bytes.length, buf.bufPos, "no trailing bytes");
    }
}
