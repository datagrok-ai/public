package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import java.sql.Timestamp;
import java.time.LocalDateTime;
import java.util.Date;

public class LocalDateTimeTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return Timestamp.valueOf((LocalDateTime) value);
    }
}
