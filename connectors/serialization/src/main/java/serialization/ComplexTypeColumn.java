package serialization;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.time.Instant;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.LocalTime;
import java.time.OffsetDateTime;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.time.temporal.Temporal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class ComplexTypeColumn extends AbstractColumn<Map<String, Object>> {
    private List<Map<String, Object>> data = new ArrayList<>();

    public ComplexTypeColumn(String name) {
        super(name);
    }

    public ComplexTypeColumn(String name, int initColumnSize) {
        super(name, initColumnSize);
    }

    @Override
    public String getType() {
        return Types.COLUMN_LIST;
    }

    @Override
    public void encode(BufferAccessor buf) {
        for (Column<?> col : convertToColumns())
            buf.writeColumn(col);
    }

    @Override
    public void add(Map<String, Object> value) {
        data.add(value != null ? value : Collections.emptyMap());
        length++;
    }

    @Override
    public void addAll(Map<String, Object>[] values) {
        for (Map<String, Object> value : values)
            add(value);
    }

    @Override
    public Map<String, Object> get(int idx) {
        return data.get(idx);
    }

    @Override
    public void set(int index, Map<String, Object> value) {
        data.set(index, value != null ? value : Collections.emptyMap());
    }

    @Override
    public long memoryInBytes() {
        return data.size() * 64L;
    }

    @Override
    public boolean isNone(int idx) {
        return idx >= data.size() || data.get(idx).isEmpty();
    }

    @Override
    public Object toArray() {
        return data.toArray(new Map[0]);
    }

    @Override
    public void empty() {
        length = 0;
        data = new ArrayList<>();
    }

    int getEncodedColumnCount() {
        return convertToColumns().length;
    }

    Column<?>[] convertToColumns() {
        if (data.isEmpty())
            return new Column<?>[0];

        // Expand rows: List<Map> values are exploded into sub-rows
        List<Map<String, Object>> expandedRows = new ArrayList<>();
        for (Map<String, Object> row : data)
            expandedRows.addAll(expandRow(row));

        // Collect all flattened key-value pairs per row, tracking column order
        Map<String, List<Object>> columns = new LinkedHashMap<>();
        int rowIndex = 0;
        for (Map<String, Object> row : expandedRows) {
            for (String key : row.keySet()) {
                if (!columns.containsKey(key)) {
                    List<Object> list = new ArrayList<>(Collections.nCopies(rowIndex, null));
                    columns.put(key, list);
                }
            }
            for (Map.Entry<String, List<Object>> entry : columns.entrySet())
                entry.getValue().add(row.get(entry.getKey()));
            rowIndex++;
        }

        List<Column<?>> result = new ArrayList<>();
        for (Map.Entry<String, List<Object>> entry : columns.entrySet()) {
            String key = entry.getKey();
            List<Object> values = entry.getValue();
            String type = detectType(values);
            String colName = name != null ? name + "." + key : key;
            Column<?> col = Column.getColumnForType(type, colName, Math.max(values.size(), initColumnSize));
            for (Object v : values)
                addValueToColumn(col, v, type);
            result.add(col);
        }
        return result.toArray(new Column<?>[0]);
    }

    /**
     * Expands a single row map into one or more sub-rows. If the map contains
     * List-of-Map values (e.g. MongoDB arrays of objects), they are "exploded":
     * scalar fields get null for the extra sub-rows, list-of-map fields fill each sub-row.
     */
    private List<Map<String, Object>> expandRow(Map<String, Object> row) {
        Map<String, Object> scalars = new LinkedHashMap<>();
        Map<String, List<Map<String, Object>>> listMaps = new LinkedHashMap<>();
        separateFields(row, "", scalars, listMaps);

        if (listMaps.isEmpty())
            return Collections.singletonList(scalars);

        int maxLen = 0;
        for (List<Map<String, Object>> list : listMaps.values())
            maxLen = Math.max(maxLen, list.size());

        List<Map<String, Object>> expanded = new ArrayList<>();
        for (int i = 0; i < maxLen; i++) {
            Map<String, Object> sub = new LinkedHashMap<>();
            // Scalar fields: value in first sub-row, null in rest
            for (Map.Entry<String, Object> entry : scalars.entrySet())
                sub.put(entry.getKey(), i == 0 ? entry.getValue() : null);
            // List-of-map fields: flatten each list item with dot-prefix
            for (Map.Entry<String, List<Map<String, Object>>> entry : listMaps.entrySet()) {
                String prefix = entry.getKey();
                List<Map<String, Object>> list = entry.getValue();
                if (i < list.size()) {
                    Map<String, Object> flat = new LinkedHashMap<>();
                    flattenMap(list.get(i), prefix, flat);
                    sub.putAll(flat);
                }
            }
            expanded.add(sub);
        }
        return expanded;
    }

    @SuppressWarnings("unchecked")
    private void separateFields(Map<String, Object> map, String prefix,
                                Map<String, Object> scalars, Map<String, List<Map<String, Object>>> listMaps) {
        for (Map.Entry<String, Object> entry : map.entrySet()) {
            String key = prefix.isEmpty() ? entry.getKey() : prefix + "." + entry.getKey();
            Object value = entry.getValue();
            if (value instanceof Map) {
                separateFields((Map<String, Object>) value, key, scalars, listMaps);
            } else if (value instanceof List && !((List<?>) value).isEmpty()
                    && ((List<?>) value).get(0) instanceof Map) {
                listMaps.put(key, (List<Map<String, Object>>) value);
            } else {
                scalars.put(key, value);
            }
        }
    }

    @SuppressWarnings("unchecked")
    private void flattenMap(Map<String, Object> map, String prefix, Map<String, Object> result) {
        for (Map.Entry<String, Object> entry : map.entrySet()) {
            String key = prefix.isEmpty() ? entry.getKey() : prefix + "." + entry.getKey();
            Object value = entry.getValue();
            if (value instanceof Map)
                flattenMap((Map<String, Object>) value, key, result);
            else
                result.put(key, value);
        }
    }

    private String detectType(List<Object> values) {
        for (Object v : values) {
            if (v == null)
                continue;
            if (v instanceof Integer || v instanceof Short || v instanceof Byte)
                return Types.INT;
            if (v instanceof Long || v instanceof BigInteger)
                return Types.BIG_INT;
            if (v instanceof Float)
                return Types.FLOAT;
            if (v instanceof Double || v instanceof BigDecimal)
                return Types.FLOAT;
            if (v instanceof Boolean)
                return Types.BOOL;
            if (v instanceof Temporal || v instanceof Date)
                return Types.DATE_TIME;
            if (v instanceof String)
                return Types.STRING;
            return Types.STRING;
        }
        return Types.STRING;
    }

    private void addValueToColumn(Column<?> col, Object value, String type) {
        if (value == null) {
            col.add(null);
            return;
        }
        switch (type) {
            case Types.INT:
                ((IntColumn) col).add(((Number) value).intValue());
                break;
            case Types.BIG_INT:
                ((BigIntColumn) col).add(value.toString());
                break;
            case Types.FLOAT:
                ((FloatColumn) col).add(((Number) value).floatValue());
                break;
            case Types.BOOL:
                ((BoolColumn) col).add((Boolean) value);
                break;
            case Types.DATE_TIME:
                ((DateTimeColumn) col).add(toEpochMicros(value));
                break;
            default:
                ((StringColumn) col).add(value.toString());
                break;
        }
    }

    private static final double MILLIS_TO_MICROS = 1000.0;

    private Double toEpochMicros(Object value) {
        long epochMillis;
        if (value instanceof Date)
            epochMillis = ((Date) value).getTime();
        else if (value instanceof LocalDate)
            epochMillis = ((LocalDate) value).atStartOfDay(ZoneId.systemDefault()).toInstant().toEpochMilli();
        else if (value instanceof LocalDateTime)
            epochMillis = ((LocalDateTime) value).atZone(ZoneId.systemDefault()).toInstant().toEpochMilli();
        else if (value instanceof Instant)
            epochMillis = ((Instant) value).toEpochMilli();
        else if (value instanceof OffsetDateTime)
            epochMillis = ((OffsetDateTime) value).toInstant().toEpochMilli();
        else if (value instanceof ZonedDateTime)
            epochMillis = ((ZonedDateTime) value).toInstant().toEpochMilli();
        else if (value instanceof LocalTime)
            epochMillis = ((LocalTime) value).atDate(LocalDate.ofEpochDay(0))
                    .atZone(ZoneId.systemDefault()).toInstant().toEpochMilli();
        else
            epochMillis = 0;
        return epochMillis * MILLIS_TO_MICROS;
    }
}
