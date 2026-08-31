package grok_connect.managers.bigint_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.resultset.ColumnMeta;
import serialization.BigIntColumn;
import serialization.Column;
import java.math.BigInteger;

public class DefaultBigIntColumnManager implements ColumnManager<Object> {
    private static final Converter<String> BIT_CONVERTER = value -> Integer.toBinaryString(new BigInteger((byte[]) value).intValue());

    @Override
    public Object convert(Object value, ColumnMeta columnMeta) {
        if (value != null && value.getClass().isArray() && value.getClass().getComponentType().equals(byte.class))
            return BIT_CONVERTER.convert(value);
        return value;
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        int precision = columnMeta.getPrecision();
        int scale = columnMeta.getScale();
        return type == java.sql.Types.BIGINT
                || typeName.equalsIgnoreCase("int8")
                || typeName.equalsIgnoreCase("serial8")
                || (type == java.sql.Types.BIT && precision > 1) || (scale > 1 && type == java.sql.Types.BIT)
                || typeName.equalsIgnoreCase("varbit");
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new BigIntColumn(name, initColumnSize);
    }
}
