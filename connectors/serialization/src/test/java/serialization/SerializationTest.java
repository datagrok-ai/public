package serialization;

import java.io.*;
import java.util.*;
import org.junit.Test;
import java.nio.file.*;
import static org.junit.Assert.*;


public class SerializationTest {
    // Path to Dart part of test.
    public static final String dart_test_path =
            "../../../ddt/bin/serialization_test/bin/serialization_test_ex.dart";

    @Test
    public void testSerialization() throws IOException, InterruptedException {
        List<Column> columns = new ArrayList<>(5);
        columns.add(new FloatColumn(new Float[] {1.0F, 2.0F, 3.0F, 4.0F}));
        columns.add(new IntColumn(new Integer[] {1, 2, 3, 4}));
        columns.add(new StringColumn(new String[] {"A", "B", "C", "D"}));
        columns.add(new BoolColumn(new Boolean[] {true, false, true, false}));
        columns.add(new DateTimeColumn(new Double[] {1099511627776000.0, 2099511627777000.0, -62135607600000000.0, -62135607600000000.0}));
        columns.add(new BigIntColumn(new String[] {"1234567890", "2345678901", "3456789012", "4567890123"}));

        columns.get(0).name = "double";
        columns.get(1).name = "int";
        columns.get(2).name = "string";
        columns.get(3).name = "bool";
        columns.get(4).name = "datetime";
        columns.get(5).name = "bigint";

        DataFrame df = new DataFrame();
        df.addColumns(columns);
        byte[] bytes = df.toByteArray();

        Files.createDirectories(Paths.get("target"));
        String path = "target/blob.bin";
        File file = new File(path);
        FileOutputStream stream = new FileOutputStream(file, false);
        stream.write(bytes, 0, bytes.length);
        stream.flush();
        stream.close();

        Process process = Runtime.getRuntime().exec("dart " + dart_test_path + " --file " + path);
        process.waitFor();

        // Print Dart test output
        String line;
        BufferedReader input = new BufferedReader(new InputStreamReader(process.getInputStream()));
        while ((line = input.readLine()) != null) {
            System.out.println(line);
        }
        input.close();

        assertEquals(0, process.exitValue());
    }
}
