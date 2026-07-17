package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.time.Instant;
import java.time.LocalDate;
import java.time.LocalTime;
import java.time.ZoneId;
import java.util.Date;

public class LocalTimeConverter implements Converter<Date> {
    private static final Logger LOGGER = LoggerFactory.getLogger(InstantTypeConverter.class);

    @Override
    public Date convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        LocalTime localTime = (LocalTime) value;
        Instant instant = localTime.atDate(LocalDate.of(1970, 1, 1))
                .atZone(ZoneId.systemDefault()).toInstant();
        // Timestamp keeps the sub-ms fraction (Date.from would truncate to ms)
        return java.sql.Timestamp.from(instant);
    }
}
