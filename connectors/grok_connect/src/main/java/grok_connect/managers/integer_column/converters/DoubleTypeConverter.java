package grok_connect.managers.integer_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DoubleTypeConverter implements Converter<Integer> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DoubleTypeConverter.class);

    @Override
    public Integer convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return ((Double) value).intValue();
    }
}
