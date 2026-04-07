package serialization;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;

@Disabled
public class SerializationTest {
    // Path to Dart part of test.
    public static final String dart_test_path =
            "../../../ddt/bin/serialization_test/bin/serialization_test_ex.dart";

    @Test
    public void testSerialization() throws IOException, InterruptedException {
        Column[] columns = new Column[6];
        columns[0] = new FloatColumn("double", new Float[] {1.0F, 2.0F, 3.0F, 4.0F});
        columns[1] = new IntColumn("int", new Integer[] {1, 2, 3, 4});
        columns[2] = new StringColumn("string", new String[] {"A", "B", "C", "D"});
        columns[3] = new BoolColumn("bool", new Boolean[] {true, false, true, false});
        columns[4] = new DateTimeColumn("datetime", new Double[] {1099511627776000.0, 2099511627777000.0, FloatColumn.None, FloatColumn.None});
        columns[5] = new BigIntColumn("bigint", new String[] {"1234567890", "2345678901", "3456789012", "4567890123"});

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

        Assertions.assertEquals(0, process.exitValue());
    }
}
