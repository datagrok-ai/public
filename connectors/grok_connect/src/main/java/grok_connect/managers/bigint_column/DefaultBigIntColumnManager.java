package grok_connect.managers.bigint_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.resultset.ColumnMeta;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.BigIntColumn;
import serialization.Column;
import java.math.BigInteger;

public class DefaultBigIntColumnManager implements ColumnManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultBigIntColumnManager.class);
    private static final Converter<String> DEFAULT_CONVERTER = Object::toString;
    private static final Converter<String> BIT_CONVERTER = value -> Integer.toBinaryString(new BigInteger((byte[]) value).intValue());

    @Override
    public String convert(Object value, ColumnMeta columnMetal) {
        LOGGER.trace("convert method was called");
        if (value == null) return "";
        if (value.getClass().isArray() && value.getClass().getComponentType().equals(byte.class))
            return BIT_CONVERTER.convert(value);
        return DEFAULT_CONVERTER.convert(value);
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
