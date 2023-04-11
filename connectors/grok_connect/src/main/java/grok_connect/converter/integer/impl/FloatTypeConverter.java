package grok_connect.converter.integer.impl;

import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FloatTypeConverter implements Converter<Integer> {
    private static final Logger LOGGER = LoggerFactory.getLogger(FloatTypeConverter.class);

    @Override
    public Integer convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return ((Float) value).intValue();
    }
}
