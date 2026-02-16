package grok_connect.managers.float_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.managers.float_column.converters.BigDecimalTypeConverter;
import grok_connect.managers.float_column.converters.DoubleTypeConverter;
import grok_connect.resultset.ColumnMeta;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.FloatColumn;
import java.math.BigDecimal;
import java.sql.Types;
import java.util.HashMap;
import java.util.Map;

public class DefaultFloatColumnManager implements ColumnManager<Float> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultFloatColumnManager.class);
    private static final Converter<Float> DEFAULT_CONVERTER = value -> (Float) value;
    private final Map<Class<?>, Converter<Float>> converterMap;

    {
        converterMap = new HashMap<>();
        converterMap.put(Double.class, new DoubleTypeConverter());
        converterMap.put(BigDecimal.class, new BigDecimalTypeConverter());
    }

    @Override
    public Float convert(Object value, ColumnMeta columnMeta) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<Float> converter = converterMap
                .getOrDefault(aClass, DEFAULT_CONVERTER);
        return converter.convert(value);
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        int scale = columnMeta.getScale();
        return type == Types.FLOAT || type == java.sql.Types.DOUBLE || type == java.sql.Types.REAL ||
                type == Types.DECIMAL ||
                typeName.equalsIgnoreCase("float8") ||
                typeName.equalsIgnoreCase("float4") ||
                typeName.equalsIgnoreCase("money") ||
                typeName.equalsIgnoreCase("binary_float") ||
                typeName.equalsIgnoreCase("binary_double") ||
                typeName.equalsIgnoreCase("numeric") ||
                typeName.equalsIgnoreCase("DECFLOAT") ||
                (typeName.equalsIgnoreCase("number") && scale > 0);
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new FloatColumn(name, initColumnSize);
    }
}
