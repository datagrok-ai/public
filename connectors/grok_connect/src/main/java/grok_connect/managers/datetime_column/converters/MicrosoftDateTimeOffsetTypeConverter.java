package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import microsoft.sql.DateTimeOffset;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Timestamp;
import java.util.Date;

public class MicrosoftDateTimeOffsetTypeConverter implements Converter<Date> {
    private static final Logger LOGGER = LoggerFactory.getLogger(MicrosoftDateTimeOffsetTypeConverter.class);

    @Override
    public Date convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        // Timestamp keeps the sub-ms fraction (Date.from would truncate to ms)
        return Timestamp.from(((DateTimeOffset) value).getOffsetDateTime().toInstant());
    }
}
