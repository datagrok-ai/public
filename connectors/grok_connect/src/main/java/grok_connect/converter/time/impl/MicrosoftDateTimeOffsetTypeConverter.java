package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import microsoft.sql.DateTimeOffset;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.Date;

public class MicrosoftDateTimeOffsetTypeConverter implements Converter<Date> {
    private static final Logger LOGGER = LoggerFactory.getLogger(MicrosoftDateTimeOffsetTypeConverter.class);

    @Override
    public Date convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return Date.from(((DateTimeOffset) value).getOffsetDateTime().toInstant());
    }
}
