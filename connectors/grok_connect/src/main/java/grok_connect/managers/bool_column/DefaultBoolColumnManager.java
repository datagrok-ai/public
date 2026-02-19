package grok_connect.managers.bool_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.resultset.ColumnMeta;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.BoolColumn;
import serialization.Column;

public class DefaultBoolColumnManager implements ColumnManager<Boolean> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultBoolColumnManager.class);
    private static final Converter<Boolean> DEFAULT_CONVERTER = value -> (Boolean) value;

    @Override
    public Boolean convert(Object value, ColumnMeta columnMeta) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        return DEFAULT_CONVERTER.convert(value);
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        int precision = columnMeta.getPrecision();
        int scale = columnMeta.getScale();
        return (type == java.sql.Types.BOOLEAN) ||
                typeName.equalsIgnoreCase("bool") || (type == java.sql.Types.BIT
                && precision == 1 && scale == 0);
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new BoolColumn(name, initColumnSize);
    }
}
