package grok_connect.converter.integer.impl;

import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ByteTypeConverter implements Converter<Integer> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ByteTypeConverter.class);

    @Override
    public Integer convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return ((Byte) value).intValue();
    }
}
