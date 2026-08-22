package serialization;

import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

// Java -> Dart golden fixtures: writes src/test/resources/d42-java/<case>.d42 with the
// production writer (DataFrame.toByteArray, best int encoders) plus an .expected.json
// sidecar in the resources/d42 shape; core/shared/ddt/test/serialization/
// java_d42_fixture_test.dart decodes every fixture with the Dart reader. Each fixture is
// also round-tripped through the Java reader here. Serialization is deterministic; the
// committed goldens are compared byte-for-byte and only rewritten with -Dd42.regenerate=true.
public class D42JavaFixtureWriterTest {
    private static final int N = IntColumn.None;
    private static final String EMOJI = "\uD83D\uDE00x";
    private static final String CJK = "\u65E5\u672C\u8A9E";
    private static final String TRICKY = "with,comma\"quote\nnewline";

    // Column spec: the column, its sidecar value encoding, the expected encoder id
    // (null = record whatever the cost model picked) and the expected per-row values.
    private static class Spec {
        final Column<?> col;
        final String enc;
        final Integer expectedId;
        final Object[] values;

        Spec(Column<?> col, String enc, Integer expectedId, Object[] values) {
            this.col = col;
            this.enc = enc;
            this.expectedId = expectedId;
            this.values = values;
        }
    }

    private static class Fixture {
        final String name;
        final List<Spec> specs = new ArrayList<>();

        Fixture(String name) {
            this.name = name;
        }

        DataFrame frame() {
            DataFrame df = new DataFrame();
            df.name = name;
            for (Spec s : specs)
                df.addColumn(s.col);
            return df;
        }
    }

    private static Object[] boxed(int[] values) {
        Object[] r = new Object[values.length];
        for (int i = 0; i < values.length; i++)
            r[i] = values[i] == N ? null : values[i];
        return r;
    }

    private static Spec ints(String name, int[] values, Integer expectedId) {
        return new Spec(new IntColumn(name, values.clone()), "int", expectedId, boxed(values));
    }

    private static Spec strings(String name, String[] values) {
        return new Spec(new StringColumn(name, values.clone()), "string", 0, values.clone());
    }

    private static Fixture intEncoders() {
        int rows = 2000;
        int[] seq = new int[rows];
        int[] blocks = new int[rows];
        int[] grouped = IntColumnEncoderSelectionTest.grouped(rows, 50);
        int[] lowcard = new int[rows];
        int[] entropy = new int[rows];
        int[] allnull = new int[rows];
        int[] walk = new int[rows];
        int[] edge = new int[rows];
        int[] rlemin = new int[rows];
        Random rnd = new Random(2026);
        int v = 0;
        for (int i = 0; i < rows; i++) {
            seq[i] = i * 3 + 7;
            blocks[i] = 1 + (i % 24) / 3;
            lowcard[i] = i % 37 == 0 ? N : rnd.nextInt(10);
            int e = rnd.nextInt();
            entropy[i] = e == N ? 0 : e;
            allnull[i] = N;
            if (i % 11 == 0)
                walk[i] = N;
            else {
                v += rnd.nextInt(5) - 2;
                walk[i] = v;
            }
            edge[i] = i % 3 == 0 ? Integer.MIN_VALUE + 1 : (i % 3 == 1 ? Integer.MAX_VALUE : N);
            // rle's none marker (min - 1) is the None sentinel itself when min == MIN_VALUE + 1.
            rlemin[i] = i % 97 == 0 ? N : Integer.MIN_VALUE + 1 + grouped[i];
        }
        Fixture fx = new Fixture("java_int_encoders");
        fx.specs.add(ints("seq", seq, 2));
        fx.specs.add(ints("blocks", blocks, 2));
        fx.specs.add(ints("grouped", grouped, 3));
        fx.specs.add(ints("lowcard", lowcard, 4));
        fx.specs.add(ints("entropy", entropy, 1));
        fx.specs.add(ints("allnull", allnull, 2));
        fx.specs.add(ints("walk", walk, null));
        fx.specs.add(ints("edge", edge, null));
        fx.specs.add(ints("rlemin", rlemin, 3));
        return fx;
    }

    private static Fixture strings() {
        int rows = 2000;
        String[] plate = new String[rows];
        String[] few = new String[rows];
        String[] unique = new String[rows];
        String[] cats = {"alpha", "beta", "gamma"};
        int[] grouped = IntColumnEncoderSelectionTest.grouped(rows, 384);
        Random rnd = new Random(5);
        for (int i = 0; i < rows; i++) {
            plate[i] = i % 997 == 0 ? null : String.format("PLT-%04d", grouped[i]);
            few[i] = i % 500 == 0 ? "" : cats[rnd.nextInt(3)];
            unique[i] = i % 3 == 0 ? "c" + i : (i % 3 == 1 ? CJK + i : EMOJI + i);
        }
        Fixture fx = new Fixture("java_strings");
        fx.specs.add(strings("plate", plate));
        fx.specs.add(strings("few", few));
        fx.specs.add(strings("unique", unique));
        return fx;
    }

    private static Fixture bigints() {
        int rows = 1000;
        BigIntColumn ids = new BigIntColumn("ids", 10);
        ids.setDowncastAllowed(true);
        BigIntColumn wide = new BigIntColumn("wide", 10);
        wide.setDowncastAllowed(true);
        Object[] idValues = new Object[rows];
        Object[] wideValues = new Object[rows];
        for (int i = 0; i < rows; i++) {
            Long id = i % 100 == 99 ? null : (long) i * 1000;
            ids.add(id);
            idValues[i] = id == null ? null : (Object) id.intValue();
            Long w = i % 7 == 0 ? null : (long) i * 3000000000L;
            wide.add(w);
            wideValues[i] = w == null ? null : w.toString();
        }
        Assertions.assertEquals(Types.INT, ids.getType());
        Assertions.assertEquals(Types.BIG_INT, wide.getType());
        Fixture fx = new Fixture("java_bigint_modes");
        fx.specs.add(new Spec(ids, "int", 3, idValues));
        fx.specs.add(new Spec(wide, "bigint", 1, wideValues));
        return fx;
    }

    private static Fixture mixed() {
        Integer[] ints = {1, null, Integer.MIN_VALUE + 1, -5, 100, 0};
        String[] strs = {"hello", "", null, CJK, EMOJI, TRICKY};
        Boolean[] bools = {true, false, true, null, false, true};
        Float[] floats = {1.5F, null, Float.NaN, -3.25F, 0.0F, 1e30F};
        Double[] doubles = {1.5, null, 1.0000000000000002, -2.0, Double.NaN, 1e308};
        Double[] dates = {0.0, -1.0, 1577836800000000.0, null, 1735108200000000.0, 1500000000000001.0};
        String[] bigs = {"1152921504606846983", null, "-42", "", "0", "999"};

        IntColumn intCol = new IntColumn("ints", ints);
        intCol.setTag("colKey", "colVal");
        Fixture fx = new Fixture("java_mixed");
        fx.specs.add(new Spec(intCol, "int", 1, ints));
        fx.specs.add(new Spec(new StringColumn("strs", strs), "string", 0, strs));
        fx.specs.add(new Spec(new BoolColumn("bools", bools), "bool", 1, bools));
        fx.specs.add(new Spec(new FloatColumn("floats", floats), "double32", 1, hex32(floats)));
        fx.specs.add(new Spec(FloatColumn.double64("doubles", doubles), "double64", 5, hex64(doubles)));
        fx.specs.add(new Spec(new DateTimeColumn("dates", dates), "datetimeMicros", 3, micros(dates)));
        fx.specs.add(new Spec(new BigIntColumn("bigs", bigs), "bigint", 1, bigs));
        return fx;
    }

    private static Object[] hex32(Float[] values) {
        Object[] r = new Object[values.length];
        for (int i = 0; i < values.length; i++)
            r[i] = values[i] == null ? null : String.format("%08x", Float.floatToRawIntBits(values[i]));
        return r;
    }

    private static Object[] hex64(Double[] values) {
        Object[] r = new Object[values.length];
        for (int i = 0; i < values.length; i++)
            r[i] = values[i] == null ? null : String.format("%016x", Double.doubleToRawLongBits(values[i]));
        return r;
    }

    private static Object[] micros(Double[] values) {
        Object[] r = new Object[values.length];
        for (int i = 0; i < values.length; i++)
            r[i] = values[i] == null ? null : (Object) values[i].longValue();
        return r;
    }

    private static Spec dates(String name, Double[] values, int expectedId) {
        return new Spec(new DateTimeColumn(name, values), "datetimeMicros", expectedId, micros(values));
    }

    // dateTime:int (id 2) for regular / jittered / date-only / None-bearing columns; microsecond columns stay on id 3.
    private static Fixture datetimes() {
        int rows = 2000;
        Fixture fx = new Fixture("java_datetimes");
        fx.specs.add(dates("regular", DateTimeColumnEncoderTest.regular(rows), 2));
        fx.specs.add(dates("jitter", DateTimeColumnEncoderTest.jitter(rows), 2));
        fx.specs.add(dates("dateonly", DateTimeColumnEncoderTest.dateOnly(rows), 2));
        fx.specs.add(dates("mixednone", DateTimeColumnEncoderTest.mixedNone(rows), 2));
        fx.specs.add(dates("subsecond", DateTimeColumnEncoderTest.subSecond(rows), 3));
        return fx;
    }

    private static Fixture zeroRow() {
        Fixture fx = new Fixture("java_zero_row");
        fx.specs.add(new Spec(new IntColumn("i", 0), "int", 1, new Object[0]));
        fx.specs.add(new Spec(new StringColumn("s", 0), "string", 0, new Object[0]));
        return fx;
    }

    private static List<Fixture> fixtures() {
        List<Fixture> list = new ArrayList<>();
        list.add(intEncoders());
        list.add(strings());
        list.add(bigints());
        list.add(mixed());
        list.add(datetimes());
        list.add(zeroRow());
        return list;
    }

    private static JsonObject sidecar(Fixture fx, Map<String, Integer> observed) {
        JsonArray columns = new JsonArray();
        for (Spec s : fx.specs) {
            JsonArray values = new JsonArray();
            for (Object v : s.values) {
                if (v == null) values.add((String) null);
                else if (v instanceof String) values.add((String) v);
                else if (v instanceof Boolean) values.add((Boolean) v);
                else values.add((Number) v);
            }
            JsonObject c = new JsonObject();
            c.addProperty("name", s.col.getName());
            c.addProperty("dgType", s.col.getType());
            c.addProperty("valueEncoding", s.enc);
            c.addProperty("encoderId", observed.get(s.col.getName()));
            if (s.col instanceof StringColumn)
                c.addProperty("indexEncoderId", IntColumnEncoderSelectionTest.indexEncoderId((StringColumn) s.col));
            c.add("values", values);
            if (!s.col.getTags().isEmpty()) {
                JsonObject tags = new JsonObject();
                for (Map.Entry<String, String> e : s.col.getTags().entrySet())
                    tags.addProperty(e.getKey(), e.getValue());
                c.add("tags", tags);
            }
            columns.add(c);
        }
        JsonObject root = new JsonObject();
        root.addProperty("frame", fx.name);
        root.addProperty("rowCount", fx.specs.get(0).col.getLength());
        root.add("columns", columns);
        return root;
    }

    private static void assertJavaRoundTrip(Fixture fx, byte[] bytes, Map<String, Integer> observed) {
        DataFrame df = DataFrame.fromByteArray(bytes);
        Assertions.assertEquals(fx.name, df.name);
        for (Spec s : fx.specs) {
            String name = s.col.getName();
            if (s.expectedId != null)
                Assertions.assertEquals(s.expectedId.intValue(), observed.get(name).intValue(), fx.name + "." + name + " encoder id");
            Column<?> col = df.getColumn(name);
            Assertions.assertEquals(s.col.getType(), col.getType(), name + " type");
            Assertions.assertEquals(s.values.length, col.getLength(), name + " length");
            for (int r = 0; r < s.values.length; r++) {
                Object e = s.values[r];
                String at = fx.name + "." + name + "[" + r + "]";
                switch (s.enc) {
                    case "int":
                        Assertions.assertEquals(e == null, col.isNone(r), at);
                        if (e != null) Assertions.assertEquals(e, col.get(r), at);
                        break;
                    case "string":
                    case "bigint": {
                        boolean none = e == null || e.toString().isEmpty();
                        Assertions.assertEquals(none, col.isNone(r), at);
                        if (!none) Assertions.assertEquals(e, String.valueOf(col.get(r)), at);
                        break;
                    }
                    case "bool":
                        Assertions.assertEquals(e != null && (Boolean) e, col.get(r), at);
                        break;
                    case "double32":
                        Assertions.assertEquals(e == null, col.isNone(r), at);
                        if (e != null) Assertions.assertEquals(e, String.format("%08x",
                                Float.floatToRawIntBits(((FloatColumn) col).get(r))), at);
                        break;
                    case "double64":
                        Assertions.assertEquals(e == null, col.isNone(r), at);
                        if (e != null) Assertions.assertEquals(e, String.format("%016x",
                                Double.doubleToRawLongBits(((FloatColumn) col).getDouble(r))), at);
                        break;
                    case "datetimeMicros":
                        Assertions.assertEquals(e == null, col.isNone(r), at);
                        if (e != null) Assertions.assertEquals(((Long) e).doubleValue(), (Double) col.get(r), 0.0, at);
                        break;
                    default:
                        Assertions.fail("unknown encoding " + s.enc);
                }
            }
            for (Map.Entry<String, String> t : s.col.getTags().entrySet())
                Assertions.assertEquals(t.getValue(), col.getTag(t.getKey()), name + " tag " + t.getKey());
        }
    }

    // Fixtures are committed under the module's own test resources (cwd = module dir under Maven).
    private static File outputDir() {
        if (!new File("src/test/java").isDirectory())
            return null;
        File dir = new File("src/test/resources/d42-java");
        dir.mkdirs();
        return dir;
    }

    private static final String REGENERATE_HINT = "; regenerate with -Dd42.regenerate=true";

    // Goldens are compared, not overwritten: a writer change must be regenerated on purpose.
    private static void checkOrWriteGolden(File file, byte[] content, boolean regenerate) throws Exception {
        if (regenerate) {
            Files.write(file.toPath(), content);
            return;
        }
        Assertions.assertTrue(file.isFile(), file + " is missing" + REGENERATE_HINT);
        byte[] onDisk = Files.readAllBytes(file.toPath());
        if (file.getName().endsWith(".json"))
            onDisk = new String(onDisk, StandardCharsets.UTF_8).replace("\r\n", "\n").getBytes(StandardCharsets.UTF_8);
        Assertions.assertArrayEquals(onDisk, content, file + " differs from the written bytes" + REGENERATE_HINT);
    }

    @Test
    public void writeFixtures() throws Exception {
        File dir = outputDir();
        boolean regenerate = Boolean.getBoolean("d42.regenerate");
        for (Fixture fx : fixtures()) {
            byte[] bytes = fx.frame().toByteArray();
            Assertions.assertArrayEquals(bytes, fx.frame().toByteArray(), fx.name + ": non-deterministic bytes");
            Assertions.assertTrue(bytes.length < 100 * 1024, fx.name + ": " + bytes.length + " bytes exceeds 100 KB");
            Map<String, Integer> observed = D42DartFixtureTest.scanEncoderIds(bytes);
            assertJavaRoundTrip(fx, bytes, observed);
            if (dir == null)
                continue;
            String json = new GsonBuilder().serializeNulls().create().toJson(sidecar(fx, observed)) + "\n";
            checkOrWriteGolden(new File(dir, fx.name + ".d42"), bytes, regenerate);
            checkOrWriteGolden(new File(dir, fx.name + ".expected.json"), json.getBytes(StandardCharsets.UTF_8), regenerate);
        }
    }
}
