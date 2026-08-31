package serialization;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import serialization.codecs.BitPacking;
import serialization.codecs.StringTokens;

// GROK-20761: Java round-trips of the level-2 codecs (string:tokens 4, int:bitPacked 5,
// int:delta 6) - the writer halves ported from ddt; cross-language proof is in
// D42DartFixtureTest via the Dart-written fixtures.
public class StringTokensCodecTest {

    private static String pad(long n, int w) {
        StringBuilder s = new StringBuilder(Long.toString(n));
        while (s.length() < w)
            s.insert(0, '0');
        return s.toString();
    }

    private static StringColumn roundTripString(String[] values, int expectedId) {
        StringColumn col = new StringColumn("s", values);
        BufferAccessor buf = new BufferAccessor();
        col.encode(buf);
        buf.bufPos = 0;
        Assertions.assertEquals(expectedId, peekId(buf));
        StringColumn copy = new StringColumn("s");
        copy.decode(buf);
        Assertions.assertEquals(values.length, copy.length);
        for (int i = 0; i < values.length; i++)
            Assertions.assertEquals(values[i] == null ? "" : values[i], copy.get(i), "row " + i);
        return copy;
    }

    private static IntColumn roundTripInt(int[] values, int expectedId) {
        IntColumn col = new IntColumn("i", values);
        Assertions.assertEquals(expectedId, IntColumn.Encoding.choose(values, values.length).id);
        BufferAccessor buf = new BufferAccessor();
        col.encode(buf);
        buf.bufPos = 0;
        Assertions.assertEquals(expectedId, peekId(buf));
        IntColumn copy = new IntColumn("i", 0);
        copy.decode(buf);
        Assertions.assertArrayEquals(values, (int[]) copy.toArray());
        return copy;
    }

    private static int peekId(BufferAccessor buf) {
        int pos = buf.bufPos;
        int id = buf.readInt32();
        buf.bufPos = pos;
        return id;
    }

    @Test
    public void bitPackingRoundTripsEveryWidth() {
        Random rnd = new Random(7);
        for (int bits = 0; bits <= 31; bits++) {
            int max = bits == 31 ? 0x7FFFFFFE : (1 << bits) - 2;
            int[] values = new int[500];
            for (int i = 0; i < values.length; i++)
                values[i] = bits == 0 ? 5 : rnd.nextInt(max + 1);
            values[3] = IntColumn.None;
            values[499] = IntColumn.None;
            BufferAccessor buf = new BufferAccessor();
            BitPacking.write(buf, values, 0, values.length);
            Assertions.assertEquals(
                    BitPacking.sizeInBytes(values.length, BitPacking.bitsFor(BitPacking.minMax(values, 0, values.length))),
                    buf.bufPos, "bits=" + bits);
            buf.bufPos = 0;
            Assertions.assertArrayEquals(values, BitPacking.read(buf), "bits=" + bits);
        }
    }

    @Test
    public void bitPackingNoneCostsABitOnlyWhenPresent() {
        Assertions.assertEquals(16, BitPacking.bitsFor(new long[]{0, 65535, 0}));
        Assertions.assertEquals(17, BitPacking.bitsFor(new long[]{0, 65535, 1}));
        Assertions.assertEquals(-1, BitPacking.bitsFor(new long[]{Integer.MIN_VALUE + 1, Integer.MAX_VALUE, 0}));
    }

    @Test
    public void intBitPackedSelectedForBoundedRandom() {
        Random rnd = new Random(7);
        int[] values = new int[5000];
        for (int i = 0; i < values.length; i++)
            values[i] = i % 97 == 0 ? IntColumn.None : rnd.nextInt(100000);
        roundTripInt(values, 5);
    }

    @Test
    public void intDeltaSelectedForSorted() {
        Random rnd = new Random(7);
        int[] values = new int[100000];
        for (int i = 0; i < values.length; i++)
            values[i] = rnd.nextInt(10000000);
        Arrays.sort(values);
        values[10] = IntColumn.None;
        roundTripInt(values, 6);
    }

    @Test
    public void tokensTemplates() {
        Assertions.assertEquals("CHEMBL-{int:6}", template("CHEMBL-000354", "CHEMBL-000355"));
        Assertions.assertEquals("CHEMBL{int}", template("CHEMBL25", "CHEMBL3"));
        Assertions.assertEquals("GO:{int:7}", template("GO:0005886", "GO:0000001"));
        Assertions.assertEquals("DSM-{int:1}:{hex:4}", template("DSM-4:555F", "DSM-4:5A5F"));
        Assertions.assertEquals("ZINC{long:12}", template("ZINC000012345678", "ZINC000012345679"));
        Assertions.assertEquals("LOT-{int:4}-{int:6}", template("LOT-2024-000123", "LOT-2025-000124"));
        Assertions.assertEquals("{int}", template("12345", "6789"));
        Assertions.assertEquals("{str}{int:2}", template("A01", "P24"));
        Assertions.assertNull(template("ABC", "DEF"));
        Assertions.assertNull(template("P12345", "Q9Y6K1", "O00123"));
        Assertions.assertNull(template("GO:0005886", "CHEMBL25", "DSM-4:555F"));
        Assertions.assertNull(template("007", "12"));
    }

    private static String template(String... cats) {
        StringTokens t = StringTokens.analyze(Arrays.asList(cats));
        return t == null ? null : t.template();
    }

    @Test
    public void tokensRoundTripsIdentifierColumns() {
        Random rnd = new Random(7);
        String[] go = new String[3000];
        for (int i = 0; i < go.length; i++)
            go[i] = i % 17 == 0 ? null : i % 23 == 0 ? "" : "GO:" + pad(rnd.nextInt(2000000), 7);
        roundTripString(go, 4);

        String[] chembl = new String[3000];
        for (int i = 0; i < chembl.length; i++)
            chembl[i] = "CHEMBL" + rnd.nextInt(5000000);
        roundTripString(chembl, 4);

        String[] hex = new String[3000];
        for (int i = 0; i < hex.length; i++) {
            String h = Integer.toHexString(rnd.nextInt(65536)).toUpperCase();
            hex[i] = "DSM-" + rnd.nextInt(9) + ":" + "0000".substring(h.length()) + h;
        }
        roundTripString(hex, 4);

        String[] zinc = new String[3000];
        for (int i = 0; i < zinc.length; i++)
            zinc[i] = "ZINC" + pad(rnd.nextInt(1000000), 6) + pad(rnd.nextInt(1000000), 6);
        roundTripString(zinc, 4);

        String[] longSplit = new String[2000];
        for (int i = 0; i < longSplit.length; i++)
            longSplit[i] = (char) (65 + rnd.nextInt(26)) + "-"
                    + pad(rnd.nextInt(1000) * 1000000000000L + rnd.nextInt(1000000), 15);
        roundTripString(longSplit, 4);

        String[] varying = new String[3000];
        for (int i = 0; i < varying.length; i++)
            varying[i] = (char) (65 + rnd.nextInt(3)) + "" + rnd.nextInt(50) + (char) (88 + rnd.nextInt(3)) + rnd.nextInt(50);
        roundTripString(varying, 4);
    }

    @Test
    public void tokensNotApplicableFallsBackToCategories() {
        Random rnd = new Random(7);
        String[] words = new String[500];
        for (int i = 0; i < words.length; i++)
            words[i] = "word" + (char) (97 + rnd.nextInt(26)) + (char) (97 + rnd.nextInt(26));
        // tokenizable ('word' + letters)? no digits -> no int part -> id 0.
        roundTripString(words, 0);

        String[] mixed = new String[500];
        for (int i = 0; i < mixed.length; i++)
            mixed[i] = i % 2 == 0 ? "GO:" + pad(rnd.nextInt(2000000), 7) : "CHEMBL" + rnd.nextInt(5000000);
        roundTripString(mixed, 0);
    }
}
