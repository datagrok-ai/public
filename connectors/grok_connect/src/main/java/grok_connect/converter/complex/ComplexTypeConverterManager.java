package grok_connect.converter.complex;

import grok_connect.GrokConnect;
import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.complex.impl.ComplexTypeConverter;
import grok_connect.resultset.ResultSetManager;
import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

public class ComplexTypeConverterManager extends AbstractConverterManager<List<Column>> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ComplexTypeConverterManager.class);
    private static final Converter<Map<String, Object>> defaultConverter = new ComplexTypeConverter();
    private static final int COLUMN_NAME_INDEX = 0;
    private final Map<String, Integer> typesMap;

    public ComplexTypeConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
        typesMap = new HashMap<>();
        typesMap.put(serialization.Types.BIG_INT, Types.BIGINT);
        typesMap.put(serialization.Types.INT, Types.INTEGER);
        typesMap.put(serialization.Types.FLOAT, Types.FLOAT);
        typesMap.put(serialization.Types.DATE_TIME, Types.DATE);
        typesMap.put(serialization.Types.BOOL, Types.BOOLEAN);
        typesMap.put(serialization.Types.STRING, Types.VARCHAR);
    }

    @Override
    public List<Column> convert(Object value, Object...args) {
        LOGGER.trace("convert method was called with args {}", Arrays.toString(args));
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        LOGGER.debug("using default converter for class {}",
                value.getClass());
        return getColumnsFromMap(convert(value), (String) args[COLUMN_NAME_INDEX], args);
    }

    protected Map<String, Object> convert(Object value) {
        return defaultConverter.convert(value);
    }

    protected List<Column> getColumnsFromMap(Map<String, Object> map, String prefix, Object... args) {
        List<Column> result = new ArrayList<>();
        String prefixFormat = "%s.%s";
        for (String key: map.keySet()) {
            Object o = map.get(key);
            if (isComplexType(o)) {
                result.addAll(getColumnsFromMap(convert(o), String.format(prefixFormat, prefix, key), args));
            } else if (o instanceof List && ((List<?>) o).get(0) instanceof Map) {
                List<?> list = (List<?>) o;
                for (Object value : list) {
                    List<Column> columnsFromMap = getColumnsFromMap(convert(value),
                            String.format(prefixFormat, prefix, key), args);
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

    protected boolean isComplexType(Object value) {
        return value instanceof Map;
    }

    protected Column getColumnForObject(Object object) {
        ResultSetManager resultSetManager = GrokConnect.providerManager.getDefaultManager();
        Column column = resultSetManager.getColumn(object);
        column.add(resultSetManager.convert(object, typesMap.get(column.getType()), "", 0, 0));
        return column;
    }
}
