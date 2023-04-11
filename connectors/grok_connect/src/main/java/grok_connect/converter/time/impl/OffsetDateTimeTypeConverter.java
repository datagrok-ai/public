package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.time.OffsetDateTime;
import java.util.Date;

public class OffsetDateTimeTypeConverter implements Converter<Date> {
    private static final Logger LOGGER = LoggerFactory.getLogger(OffsetDateTimeTypeConverter.class);

    @Override
    public Date convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return Date.from(((OffsetDateTime) value)
                .toInstant());
    }
}
