package grok_connect.converter.complex;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.complex.impl.ComplexTypeConverter;
import grok_connect.resultset.ResultSetManager;
import grok_connect.type.TypeChecker;
import microsoft.sql.DateTimeOffset;
import oracle.sql.DATE;
import oracle.sql.TIMESTAMP;
import oracle.sql.TIMESTAMPTZ;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.*;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.time.temporal.Temporal;
import java.util.*;

public class ComplexTypeConverterManager extends AbstractConverterManager<List<Column>> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ComplexTypeConverterManager.class);
    private static final Converter<Map<String, Object>> defaultConverter = new ComplexTypeConverter();
    private static final int COLUMN_NAME_INDEX = 0;
    private static final int RESULT_SET_MANAGER_INDEX = 1;
    public ComplexTypeConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
    }

    @Override
    public List<Column> convert(Object value, Object...args) {
        LOGGER.trace("convert method was called with args");
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
                Column columnForObject = getColumnForObject(o, (ResultSetManager) args[RESULT_SET_MANAGER_INDEX]);
                columnForObject.name = String.format(prefixFormat, prefix, key);
                result.add(columnForObject);
            }
        }
        return result;
    }

    protected boolean isComplexType(Object value) {
        return value instanceof Map;
    }

    protected Column getColumnForObject(Object object, ResultSetManager resultSetManager) {
        Column column = new StringColumn();
        if (object instanceof Byte || object instanceof Short || object instanceof Integer) {
            column = new IntColumn();
            column.add(resultSetManager.convert(object, 0, "int", 0, 0));
        } else if (object instanceof Long || object instanceof BigInteger) {
            column = new BigIntColumn();
            column.add(resultSetManager.convert(object, -5, "bigint", 0, 0));
        } else if (object instanceof Float || object instanceof Double || object instanceof BigDecimal) {
            column = new FloatColumn();
            column.add(resultSetManager.convert(object, 0, "float8", 0, 0));
        } else if (object instanceof Boolean) {
            column = new BoolColumn();
            column.add(object);
        } else if (object instanceof Temporal || object instanceof Date
                || object instanceof DateTimeOffset || object instanceof DATE
                || object instanceof TIMESTAMP || object instanceof TIMESTAMPTZ) {
            column = new DateTimeColumn();
            column.add(resultSetManager.convert(object, 91, "", 0, 0));
        } else {
            column.add(object.toString());
        }
        return column;
    }
}
