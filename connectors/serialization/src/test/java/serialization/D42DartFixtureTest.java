package serialization;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.io.File;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

// GROK-20345 (WO-2): cross-language proof - Dart-written d42 blobs decode in the
// Java WO-1 reader value-for-value. Fixtures + sidecars under resources/d42 are
// byte-identical copies of the Dart generator output (core/shared/ddt
// test/serialization/d42_fixture_generator_test.dart; regeneration in the README).
//
// Each column's observed on-wire encoder id is asserted against the sidecar so a
// ddt cost-model drift fails loudly; float64 is compared bit-exact, float32 at
// float precision (raw bits). testIdInventoryCovered pins the WO-1 decoder set.
public class D42DartFixtureTest {

    private static File fixtureDir() {
        try {
            URL url = D42DartFixtureTest.class.getResource("/d42");
            if (url == null)
                throw new IllegalStateException("d42 fixtures not found on the test classpath");
            return new File(url.toURI());
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    static Stream<String> fixtures() {
        File[] files = fixtureDir().listFiles((d, n) -> n.endsWith(".d42"));
        Assertions.assertNotNull(files, "no d42 fixtures");
        List<String> names = new ArrayList<>();
        for (File f : files)
            names.add(f.getName().substring(0, f.getName().length() - ".d42".length()));
        Collections.sort(names);
        Assertions.assertFalse(names.isEmpty(), "no d42 fixtures");
        return names.stream();
    }

    private static JsonObject readSidecar(File dir, String base) throws Exception {
        String json = new String(Files.readAllBytes(new File(dir, base + ".expected.json").toPath()),
                StandardCharsets.UTF_8);
        return JsonParser.parseString(json).getAsJsonObject();
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("fixtures")
    public void fixtureDecodes(String base) throws Exception {
        File dir = fixtureDir();
        byte[] bytes = Files.readAllBytes(new File(dir, base + ".d42").toPath());
        JsonObject sidecar = readSidecar(dir, base);

        DataFrame df = DataFrame.fromByteArray(bytes);
        Map<String, Integer> observed = scanEncoderIds(bytes);

        Assertions.assertEquals(sidecar.get("frame").getAsString(), df.name, base + ": frame name");
        int rowCount = sidecar.get("rowCount").getAsInt();
        Assertions.assertEquals(rowCount, df.rowCount.intValue(), base + ": rowCount");

        for (JsonElement ce : sidecar.getAsJsonArray("columns")) {
            JsonObject c = ce.getAsJsonObject();
            String name = c.get("name").getAsString();
            String enc = c.get("valueEncoding").getAsString();
            int expectedId = c.get("encoderId").getAsInt();

            Column<?> col = df.getColumn(name);
            Assertions.assertNotNull(col, base + ": missing column " + name);
            Assertions.assertEquals(expectedId, observed.get(name).intValue(),
                    base + "." + name + ": observed encoder id != sidecar");

            JsonArray values = c.getAsJsonArray("values");
            Assertions.assertEquals(rowCount, values.size(), base + "." + name + ": value count");
            for (int r = 0; r < rowCount; r++) {
                JsonElement v = values.get(r);
                boolean isNull = v.isJsonNull();
                String at = base + "." + name + "[" + r + "]";
                switch (enc) {
                    case "int":
                        assertNone(col, r, isNull, at);
                        if (!isNull) Assertions.assertEquals(Integer.valueOf(v.getAsInt()),
                                ((IntColumn) col).get(r), at);
                        break;
                    case "double32": {
                        FloatColumn f = (FloatColumn) col;
                        if (isNull) Assertions.assertTrue(f.isNone(r), at + " expected none");
                        else {
                            Assertions.assertFalse(f.isNone(r), at + " unexpected none");
                            int expBits = Integer.parseUnsignedInt(v.getAsString(), 16);
                            Assertions.assertEquals(expBits, Float.floatToRawIntBits(f.get(r)), at + " float32 bits");
                        }
                        break;
                    }
                    case "double64": {
                        FloatColumn f = (FloatColumn) col;
                        if (isNull) Assertions.assertTrue(f.isNone(r), at + " expected none");
                        else {
                            Assertions.assertFalse(f.isNone(r), at + " unexpected none");
                            long expBits = Long.parseUnsignedLong(v.getAsString(), 16);
                            Assertions.assertEquals(expBits, Double.doubleToRawLongBits(f.getDouble(r)), at + " float64 bits");
                        }
                        break;
                    }
                    case "bool":
                        // bool is never none; null decodes to false.
                        Assertions.assertEquals(!isNull && v.getAsBoolean(), ((BoolColumn) col).get(r), at);
                        break;
                    case "string": {
                        boolean none = isNull || v.getAsString().isEmpty();
                        Assertions.assertEquals(none, col.isNone(r), at + " none");
                        if (!none) Assertions.assertEquals(v.getAsString(), ((StringColumn) col).get(r), at);
                        break;
                    }
                    case "bigint": {
                        boolean none = isNull || v.getAsString().isEmpty();
                        Assertions.assertEquals(none, col.isNone(r), at + " none");
                        if (!none) Assertions.assertEquals(v.getAsString(), ((BigIntColumn) col).get(r), at);
                        break;
                    }
                    case "datetimeMicros":
                        assertNone(col, r, isNull, at);
                        if (!isNull) Assertions.assertEquals(v.getAsDouble(),
                                ((DateTimeColumn) col).get(r), 0.0, at);
                        break;
                    default:
                        Assertions.fail(base + ": unknown valueEncoding " + enc);
                }
            }
        }
    }

    private static void assertNone(Column<?> col, int r, boolean expectedNone, String at) {
        Assertions.assertEquals(expectedNone, col.isNone(r), at + " none");
    }

    @Test
    public void testIdInventoryCovered() throws Exception {
        File dir = fixtureDir();
        File[] sidecars = dir.listFiles((d, n) -> n.endsWith(".expected.json"));
        Assertions.assertNotNull(sidecars);

        Map<String, Set<Integer>> seen = new HashMap<>();
        for (File f : sidecars) {
            String base = f.getName().substring(0, f.getName().length() - ".expected.json".length());
            JsonObject sc = readSidecar(dir, base);
            for (JsonElement ce : sc.getAsJsonArray("columns")) {
                JsonObject c = ce.getAsJsonObject();
                seen.computeIfAbsent(c.get("dgType").getAsString(), k -> new HashSet<>())
                        .add(c.get("encoderId").getAsInt());
            }
        }

        Map<String, int[]> required = new LinkedHashMap<>();
        required.put("int", new int[]{1, 2, 3, 4});
        required.put("string", new int[]{0});
        required.put("bool", new int[]{1});
        required.put("datetime", new int[]{3});
        required.put("double", new int[]{1, 5});
        required.put("bigint", new int[]{1});

        for (Map.Entry<String, int[]> e : required.entrySet())
            for (int id : e.getValue())
                Assertions.assertTrue(
                        seen.getOrDefault(e.getKey(), Collections.emptySet()).contains(id),
                        "no committed fixture exercises " + e.getKey() + " encoder id " + id
                                + " (seen: " + seen.get(e.getKey()) + ")");
    }

    // Re-parses the DataFrame region (single-table blob starts at offset 0) to
    // recover each column's on-wire encoder id - mirrors the Dart scan in the
    // generator. Peeks the id (the Int32 after the column header), then decodes to
    // advance to the next column.
    static Map<String, Integer> scanEncoderIds(byte[] bytes) {
        BufferAccessor buf = new BufferAccessor(bytes);
        buf.bufPos = 0;
        buf.readTypeCode(BufferAccessor.TYPE_DATA_FRAME);
        buf.readInt64(); // rowCount
        int colCount = (int) buf.readInt64();
        buf.readString(); // name
        buf.readStringMap(); // tags

        Map<String, Integer> ids = new LinkedHashMap<>();
        for (int i = 0; i < colCount; i++) {
            buf.readTypeCode(BufferAccessor.TYPE_COLUMN);
            String name = buf.readString();
            String type = buf.readString();
            buf.readStringMap();
            int idPos = buf.bufPos;
            int encoderId = buf.readInt32();
            buf.bufPos = idPos; // rewind so decode reads the id itself
            Column<?> col = Column.getColumnForType(type, name, 0);
            col.decode(buf);
            ids.put(name, encoderId);
        }
        return ids;
    }
}
