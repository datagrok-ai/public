package grok_connect.column.bigint;

import grok_connect.column.ColumnManager;
import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.BigIntColumn;
import serialization.Column;
import java.math.BigInteger;

public class DefaultBigIntColumnManager implements ColumnManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultBigIntColumnManager.class);
    private static final Converter<String> DEFAULT_CONVERTER = Object::toString;

    @Override
    public String convert(Object value, Object... args) {
        LOGGER.trace("convert method was called");
        LOGGER.trace("using default converter");
        return value == null ? "" : DEFAULT_CONVERTER.convert(value);
    }

    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return type == java.sql.Types.BIGINT
                || typeName.equalsIgnoreCase("int8")
                || typeName.equalsIgnoreCase("serial8")
                || (type == java.sql.Types.BIT && precision > 1) || (scale > 1 && type == java.sql.Types.BIT)
                || typeName.equalsIgnoreCase("varbit");
    }

    @Override
    public boolean isApplicable(Object o) {
        return o instanceof Long || o instanceof BigInteger;
    }

    @Override
    public Column getColumn() {
        return new BigIntColumn();
    }

    @Override
    public Column getColumnWithInitSize(int size) {
        return new BigIntColumn(new String[size]);
    }
}
