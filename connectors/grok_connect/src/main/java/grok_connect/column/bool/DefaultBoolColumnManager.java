package grok_connect.column.bool;

import grok_connect.column.ColumnManager;
import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.BoolColumn;
import serialization.Column;

public class DefaultBoolColumnManager implements ColumnManager<Boolean> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultBoolColumnManager.class);
    private static final Converter<Boolean> DEFAULT_CONVERTER = value -> (Boolean) value;

    @Override
    public Boolean convert(Object value, Object... args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        return DEFAULT_CONVERTER.convert(value);
    }

    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.BOOLEAN) ||
                typeName.equalsIgnoreCase("bool") || (type == java.sql.Types.BIT
                && precision == 1 && scale == 0);
    }

    @Override
    public boolean isApplicable(Object o) {
        return o instanceof Boolean;
    }

    @Override
    public Column getColumn() {
        return new BoolColumn();
    }

    @Override
    public Column getColumnWithInitSize(int size) {
        return new BoolColumn(new Boolean[size]);
    }
}
