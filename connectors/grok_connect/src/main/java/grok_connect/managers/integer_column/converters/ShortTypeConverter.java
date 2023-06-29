package grok_connect.managers.integer_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ShortTypeConverter implements Converter<Integer> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ShortTypeConverter.class);

    @Override
    public Integer convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return ((Short)value).intValue();
    }
}
