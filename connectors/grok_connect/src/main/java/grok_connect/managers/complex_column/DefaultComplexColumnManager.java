package grok_connect.managers.complex_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.managers.complex_column.converters.ComplexTypeConverter;
import grok_connect.resultset.ColumnMeta;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import microsoft.sql.DateTimeOffset;
import oracle.sql.DATE;
import oracle.sql.TIMESTAMP;
import oracle.sql.TIMESTAMPTZ;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.ComplexTypeColumn;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.sql.Types;
import java.time.temporal.Temporal;
import java.util.*;

public class DefaultComplexColumnManager implements ColumnManager<List<Column>> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultComplexColumnManager.class);
    private static final Converter<Map<String, Object>> DEFAULT_CONVERTER = new ComplexTypeConverter();
    private final Map<String, Integer> typesMap;
    private final Map<Class<?>, String> columnsMap;
    private ColumnMeta columnMeta;
    ResultSetManager resultSetManager;

    {
        typesMap = new HashMap<>();
        typesMap.put(serialization.Types.BIG_INT, Types.BIGINT);
        typesMap.put(serialization.Types.INT, Types.INTEGER);
        typesMap.put(serialization.Types.FLOAT, Types.FLOAT);
        typesMap.put(serialization.Types.DATE_TIME, Types.DATE);
        typesMap.put(serialization.Types.BOOL, Types.BOOLEAN);
        typesMap.put(serialization.Types.STRING, Types.VARCHAR);

        columnsMap = new HashMap<>();
        columnsMap.put(Long.class, serialization.Types.BIG_INT);
        columnsMap.put(BigInteger.class, serialization.Types.BIG_INT);
        columnsMap.put(Boolean.class, serialization.Types.BOOL);
        columnsMap.put(Map.class, serialization.Types.COLUMN_LIST);
        columnsMap.put(Temporal.class, serialization.Types.DATE_TIME);
        columnsMap.put(Date.class, serialization.Types.DATE_TIME);
        columnsMap.put(java.sql.Date.class, serialization.Types.DATE_TIME);
        columnsMap.put(DateTimeOffset.class, serialization.Types.DATE_TIME);
        columnsMap.put(DATE.class, serialization.Types.DATE_TIME);
        columnsMap.put(TIMESTAMP.class, serialization.Types.DATE_TIME);
        columnsMap.put(TIMESTAMPTZ.class, serialization.Types.DATE_TIME);
        columnsMap.put(Float.class, serialization.Types.FLOAT);
        columnsMap.put(Double.class, serialization.Types.FLOAT);
        columnsMap.put(BigDecimal.class, serialization.Types.FLOAT);
        columnsMap.put(Byte.class, serialization.Types.INT);
        columnsMap.put(Short.class, serialization.Types.INT);
        columnsMap.put(Integer.class, serialization.Types.INT);
        columnsMap.put(String.class, serialization.Types.STRING);
    }

    @Override
    public List<Column> convert(Object value, ColumnMeta columnMeta) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null, return empty list");
            return new ArrayList<>();
        }
        LOGGER.debug("using default converter for class {}",
                value.getClass());
        this.columnMeta = columnMeta;
        return getColumnsFromMap(convert(value), columnMeta.getColumnLabel());
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return false;
    }

    @Override
    public Column getColumn() {
        return new ComplexTypeColumn();
    }

    @Override
    public Column getColumn(int initColumnSize) {
        return new ComplexTypeColumn(initColumnSize);
    }

    protected Map<String, Object> convert(Object value) {
        return DEFAULT_CONVERTER.convert(value);
    }

    protected List<Column> getColumnsFromMap(Map<String, Object> map, String prefix) {
        List<Column> result = new ArrayList<>();
        String prefixFormat = "%s.%s";
        for (String key: map.keySet()) {
            Object o = map.get(key);
            if (isApplicable(o)) {
                result.addAll(getColumnsFromMap(convert(o), String.format(prefixFormat, prefix, key)));
            } else if (o instanceof List && ((List<?>) o).get(0) instanceof Map) {
                List<?> list = (List<?>) o;
                for (Object value : list) {
                    List<Column> columnsFromMap = getColumnsFromMap(convert(value),
                            String.format(prefixFormat, prefix, key));
                    for (Column column: columnsFromMap) {
                        Optional<Column> foundColumn= result.stream().filter(column1 -> column1.name.equals(column.name)).findFirst();
                        if (foundColumn.isPresent()) {
                            foundColumn.get().add(column.get(0));
                        } else {
                            result.add(column);
                        }
                    }
                }
            } else {
                Column columnForObject = getColumnForObject(o);
                columnForObject.name = String.format(prefixFormat, prefix, key);
                result.add(columnForObject);
            }
        }
        return result;
    }

    protected Column getColumnForObject(Object object) {
        if (resultSetManager == null)
            resultSetManager = DefaultResultSetManager.getDefaultManager();
        Column column = Column.getColumnForType(columnsMap.getOrDefault(object.getClass(), serialization.Types.STRING), columnMeta.getColumnSize());
        ColumnMeta columnMeta = new ColumnMeta(typesMap.get(column.getType()), object.getClass().getSimpleName(), 0, 0,
                "", this.columnMeta.getColumnSize());
        column.add(resultSetManager.getApplicableColumnManager(columnMeta).convert(object, columnMeta));
        return column;
    }

    private boolean isApplicable(Object o) {
        return o instanceof Map;
    }
}
