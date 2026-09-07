package serialization;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

// Column offsets and encoder ids exposed for the GrokConnect "COLUMN SIZES" debug line.
public class TablesBlobColumnOffsetsTest {

    @Test
    public void offsetsAndEncoderIdsCoverEveryColumn() {
        DataFrame df = new DataFrame();
        df.addColumn(new IntColumn("ints", new Integer[] {1, 2, null, 4}));
        df.addColumn(new StringColumn("strs", new String[] {"a", "b", "a", null}));
        df.addColumn(new FloatColumn("floats", new Float[] {1.5F, null, 2.5F, 3.5F}));

        TablesBlob blob = df.toBlob();
        byte[] bytes = blob.toByteArray();
        Assertions.assertArrayEquals(bytes, df.toByteArray());

        int[] offsets = blob.getColumnsOffsets(0);
        Assertions.assertEquals(3, offsets.length);
        Assertions.assertTrue(offsets[0] > 0);
        Assertions.assertTrue(offsets[0] < offsets[1] && offsets[1] < offsets[2]);
        Assertions.assertTrue(offsets[2] < blob.getLength() && blob.getLength() < bytes.length);

        int[] expectedIds = {4, 0, 1};
        BufferAccessor reader = new BufferAccessor(bytes);
        int total = 0;
        for (int i = 0; i < offsets.length; i++) {
            int end = i + 1 < offsets.length ? offsets[i + 1] : blob.getLength();
            total += end - offsets[i];
            reader.bufPos = offsets[i];
            Assertions.assertEquals(expectedIds[i], reader.readColumnEncoderId(), "encoder id of column " + i);
        }
        Assertions.assertEquals(blob.getLength() - offsets[0], total);

        reader.bufPos = offsets[1];
        Assertions.assertEquals("strs", ((StringColumn) reader.readColumn()).getName());
        Assertions.assertEquals(offsets[2], reader.bufPos);
    }
}
