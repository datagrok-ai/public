package serialization;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class StringColumn extends AbstractColumn<String> {
    private static final String TYPE = Types.STRING;

    private String[] data;
    private Integer[] idxs;
    private List<String> categories;

    public StringColumn(String name) {
        super(name);
        data = new String[initColumnSize];
    }

    public StringColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
        data = new String[initColumnSize];
    }

    public StringColumn(String name, String[] values) {
        super(name);
        data = new String[initColumnSize];
        addAll(values);
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public void empty() {
        length = 0;
        data = new String[initColumnSize];
        categorize();
    }

    @Override
    public void encode(BufferAccessor buf) {
        categorize();
        buf.writeInt32(0);
        buf.writeStringList(categories.toArray(new String[0]));
        IntColumn col = new IntColumn("");
        col.addAll(idxs);
        col.encode(buf);
    }

    @Override
    public void add(String value) {
        ensureSpace(1);
        data[length++] = value;
    }

    @Override
    public void addAll(String[] values) {
        ensureSpace(values.length);
        for (String value : values)
            data[length++] = value;
    }

    @Override
    public String get(int idx) {
        return data[idx];
    }

    @Override
    public void set(int index, String value) {
        data[index] = value;
    }

    @Override
    public long memoryInBytes() {
        long size = 0;
        for (String s : data)
            if (s != null)
                size += s.length();
        return size * 2 + (long) data.length * 16;
    }

    @Override
    public boolean isNone(int idx) {
        return data[idx] == null || data[idx].equals("");
    }

    @Override
    public Object toArray() {
        return data;
    }

    private void ensureSpace(int extraLength) {
        if (length + extraLength > data.length) {
            String[] newData = new String[data.length * 2 + Math.max(0, length + extraLength - data.length * 2)];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
        }
    }

    private void categorize() {
        Map<String, Integer> categoryMap = new HashMap<>();
        int[] tempIdxs = new int[length];

        for (int n = 0; n < length; n++) {
            String value = data[n] == null ? "" : data[n];
            data[n] = value;
            tempIdxs[n] = categoryMap.computeIfAbsent(value, k -> categoryMap.size());
        }

        categories = new ArrayList<>(categoryMap.keySet());
        categories.sort((a, b) -> {
            if (a.isEmpty() && b.isEmpty()) return 0;
            if (a.isEmpty()) return 1;
            if (b.isEmpty()) return -1;
            return a.compareTo(b);
        });

        int[] remap = new int[categoryMap.size()];
        for (int i = 0; i < categories.size(); i++)
            remap[categoryMap.get(categories.get(i))] = i;

        idxs = new Integer[length];
        for (int n = 0; n < length; n++)
            idxs[n] = remap[tempIdxs[n]];
    }
}
