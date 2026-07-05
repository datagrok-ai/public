package serialization;

import java.util.ArrayList;
import java.util.List;

// Read-capable minimal port of DataFrameColumn (dataFrame:raw, encoder id 1).
// Storage is List<DataFrame>; a null entry (null nested frame) is None. Ports
// data_frame_column_encoders.dart:13-21 (recursive readDataFrame per cell). The
// Java writer never produced these and the mutation transport does not use them -
// the reader must not choke on them.
public class DataFrameColumn extends AbstractColumn<DataFrame> {
    private static final String TYPE = Types.DATA_FRAME;

    private List<DataFrame> data = new ArrayList<>();

    public DataFrameColumn(String name) {
        super(name);
    }

    public DataFrameColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        data = new ArrayList<>();
    }

    @Override
    public void encode(BufferAccessor buf) {
        throw new UnsupportedOperationException("DataFrameColumn encode is not supported");
    }

    @Override
    public void decode(BufferAccessor buf) {
        int id = buf.readInt32();
        if (id != 1)
            throw new RuntimeException("decoding " + name + ": dataframe encoder " + id + " not found");
        int len = buf.readInt32();
        data = new ArrayList<>(Math.max(len, 0));
        for (int i = 0; i < len; i++)
            data.add(buf.readDataFrame(null));
        length = len;
    }

    @Override
    public void add(DataFrame value) {
        data.add(value);
        length++;
    }

    @Override
    public void addAll(DataFrame[] values) {
        for (DataFrame value : values)
            add(value);
    }

    @Override
    public DataFrame get(int idx) {
        return data.get(idx);
    }

    @Override
    public void set(int index, DataFrame value) {
        data.set(index, value);
    }

    @Override
    public long memoryInBytes() {
        long sum = 0;
        for (DataFrame df : data)
            sum += (df == null) ? 0 : df.memoryInBytes();
        return sum;
    }

    @Override
    public boolean isNone(int idx) {
        return data.get(idx) == null;
    }

    @Override
    public Object toArray() {
        return data.toArray(new DataFrame[0]);
    }
}
