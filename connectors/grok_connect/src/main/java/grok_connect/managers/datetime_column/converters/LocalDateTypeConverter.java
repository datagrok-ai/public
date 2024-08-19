package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.time.LocalDate;
import java.time.ZoneId;
import java.util.Date;

public class LocalDateTypeConverter implements Converter<Date> {
    private static final Logger LOGGER = LoggerFactory.getLogger(LocalDateTypeConverter.class);

    @Override
    public Date convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        return Date.from(((LocalDate) value).atStartOfDay(ZoneId.systemDefault())
                .toInstant());
    }
}
