package serialization;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;


// String column.
public class StringColumn extends Column<String> {
    private static final String TYPE = Types.STRING;

    private String[] data;
    private Integer[] idxs;
    private List<String> categories;

    public StringColumn() {
        data = new String[initColumnSize];
    }

    public StringColumn(String name) {
        this();
        this.name = name;
    }

    public StringColumn(int initColumnSize) {
        this.initColumnSize = initColumnSize;
        data = new String[initColumnSize];
    }

    public StringColumn(String[] values) {
        data = new String[initColumnSize];
        addAll(values);
    }

    public String getType() {
        return TYPE;
    }

    public void empty() {
        length = 0;
        data = new String[initColumnSize];
        categorize();
    }

    public void encode(BufferAccessor buf) {
        categorize();
        buf.writeInt32(0);
        buf.writeStringList(categories.toArray(new String[0]));
        IntColumn col = new IntColumn();
        col.addAll(idxs);
        col.encode(buf);
    }

    public void add(String value) {
        ensureSpace(1);
        data[length++] = value;
    }

    public void addAll(String[] values) {
        ensureSpace(values.length);
        for (int n = 0; n < values.length; n++)
            data[length++] = values[n];
    }

    public Object get(int idx) {
        return data[idx];
    }

    @Override
    public void set(int index, String value) {
        data[index] = value;
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            String[] newData = new String[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    public int comparer(Integer o1, Integer o2) {
        if (data[o1].equals("") && data[o2].equals("")) return 0;
        if (data[o1].equals("")) return 1;
        if (data[o2].equals("")) return -1;
        return data[o1].compareTo(data[o2]);
    }

    @Override
    public long memoryInBytes() {
        long size = 0;
        for (int n = 0; n < data.length; n++)
            if (data[n] != null)
                size += data[n].length();
        return size * 2 + data.length * 16;
    }

    public boolean isNone(int idx) {
        return data[idx] == null || data[idx].equals("");
    }

    private void categorize() {
        categories = new ArrayList<>();

        idxs = new Integer[length];
        Integer[] order = new Integer[length];
        for (int n = 0; n < length; n++) {
            data[n] = data[n] == null ? "" : data[n]; 
            idxs[n] = n;
            order[n] = n;
        }

        Arrays.sort(order, new Comparator<Integer>() {
            public int compare(Integer o1, Integer o2) {
                return comparer(o1, o2);
            }
        });

        for (int i = 0; i < length; i++) {
            boolean newCat = (i == 0) || (comparer(order[i], order[i - 1]) != 0);
            if (newCat) {
                String cat = data[order[i]];
                categories.add(cat);
            }
            idxs[order[i]] = categories.size() - 1;
        }
    }

    public String[] getData() {
        return data;
    }
}
