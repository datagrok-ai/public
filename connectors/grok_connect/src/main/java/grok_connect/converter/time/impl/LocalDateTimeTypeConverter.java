package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import java.sql.Timestamp;
import java.time.LocalDateTime;
import java.util.Date;

public class LocalDateTimeTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return Timestamp.valueOf((LocalDateTime) value);
    }
}
