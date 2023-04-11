package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Timestamp;
import java.time.LocalDateTime;
import java.util.Date;

public class LocalDateTimeTypeConverter implements Converter<Date> {
    private static final Logger LOGGER = LoggerFactory.getLogger(LocalDateTimeTypeConverter.class);

    @Override
    public Date convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return Timestamp.valueOf((LocalDateTime) value);
    }
}
