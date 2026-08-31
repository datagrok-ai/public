package grok_connect.managers.float_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.managers.float_column.converters.BigDecimalTypeConverter;
import grok_connect.managers.float_column.converters.DoubleTypeConverter;
import grok_connect.resultset.ColumnMeta;
import serialization.Column;
import serialization.FloatColumn;
import java.math.BigDecimal;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Types;
import java.util.HashMap;
import java.util.Map;

public class DefaultFloatColumnManager implements ColumnManager<Float> {
    private static final Converter<Float> DEFAULT_CONVERTER = value -> (Float) value;
    private final Map<Class<?>, Converter<Float>> converterMap;

    {
        converterMap = new HashMap<>();
        converterMap.put(Double.class, new DoubleTypeConverter());
        converterMap.put(BigDecimal.class, new BigDecimalTypeConverter());
    }

    @Override
    public Float convert(Object value, ColumnMeta columnMeta) {
        if (value == null)
            return null;
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
    public boolean canReadFast(ColumnMeta columnMeta) {
        return isDouble(columnMeta) || columnMeta.getType() == Types.REAL
                || columnMeta.getTypeName().equalsIgnoreCase("float4");
    }

    // float8 keeps getObject's double -> float narrowing: getFloat would round the text once and
    // can differ from the narrowed double by an ulp.
    @Override
    public void readFast(ResultSet resultSet, int index, Column<?> column, ColumnMeta columnMeta) throws SQLException {
        float value = isDouble(columnMeta) ? DoubleTypeConverter.narrow(resultSet.getDouble(index)) : resultSet.getFloat(index);
        if (resultSet.wasNull())
            ((FloatColumn) column).add(null);
        else
            ((FloatColumn) column).add(value);
    }

    private static boolean isDouble(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        return type == Types.DOUBLE || type == Types.FLOAT || columnMeta.getTypeName().equalsIgnoreCase("float8");
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new FloatColumn(name, initColumnSize);
    }
}
