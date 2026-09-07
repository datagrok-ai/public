package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import java.time.Instant;
import java.time.LocalDate;
import java.time.LocalTime;
import java.time.ZoneId;
import java.util.Date;

public class LocalTimeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        LocalTime localTime = (LocalTime) value;
        Instant instant = localTime.atDate(LocalDate.of(1970, 1, 1))
                .atZone(ZoneId.systemDefault()).toInstant();
        // Timestamp keeps the sub-ms fraction (Date.from would truncate to ms)
        return java.sql.Timestamp.from(instant);
    }
}
