package grok_connect.managers.complex_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.managers.complex_column.converters.ComplexTypeConverter;
import grok_connect.resultset.ColumnMeta;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.ComplexTypeColumn;
import java.util.Collections;
import java.util.Map;

public class DefaultComplexColumnManager implements ColumnManager<Map<String, Object>> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultComplexColumnManager.class);
    private static final Converter<Map<String, Object>> DEFAULT_CONVERTER = new ComplexTypeConverter();

    @Override
    public Map<String, Object> convert(Object value, ColumnMeta columnMeta) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null, return empty map");
            return Collections.emptyMap();
        }
        LOGGER.debug("using default converter for class {}",
                value.getClass());
        return DEFAULT_CONVERTER.convert(value);
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return false;
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new ComplexTypeColumn(name, initColumnSize);
    }
}
