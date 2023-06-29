package grok_connect.managers.complex_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.managers.complex_column.converters.ComplexTypeConverter;
import grok_connect.resultset.ColumnMeta;
import grok_connect.resultset.DefaultResultSetManager;
import grok_connect.resultset.ResultSetManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.ComplexTypeColumn;
import java.sql.Types;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

public class DefaultComplexColumnManager implements ColumnManager<List<Column>> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultComplexColumnManager.class);
    private static final Converter<Map<String, Object>> DEFAULT_CONVERTER = new ComplexTypeConverter();
    private final Map<String, Integer> typesMap;

    {
        typesMap = new HashMap<>();
        typesMap.put(serialization.Types.BIG_INT, Types.BIGINT);
        typesMap.put(serialization.Types.INT, Types.INTEGER);
        typesMap.put(serialization.Types.FLOAT, Types.FLOAT);
        typesMap.put(serialization.Types.DATE_TIME, Types.DATE);
        typesMap.put(serialization.Types.BOOL, Types.BOOLEAN);
        typesMap.put(serialization.Types.STRING, Types.VARCHAR);
    }

    @Override
    public List<Column> convert(Object value, String columnLabel) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null, return empty list");
            return new ArrayList<>();
        }
        LOGGER.debug("using default converter for class {}",
                value.getClass());
        return getColumnsFromMap(convert(value), columnLabel);
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return false;
    }

    @Override
    public boolean isApplicable(Object o) {
        return o instanceof Map;
    }

    @Override
    public Column getColumn() {
        return new ComplexTypeColumn();
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
        ResultSetManager resultSetManager = DefaultResultSetManager.getDefaultManager();
        Column column = resultSetManager.getColumn(object);
        ColumnMeta columnMeta = new ColumnMeta(typesMap.get(column.getType()), object.getClass().getSimpleName(), 0, 0,
                "");
        column.add(resultSetManager.convert(object, columnMeta));
        return column;
    }
}
