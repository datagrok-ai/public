package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import java.time.LocalDate;
import java.time.ZoneId;
import java.util.Date;

public class LocalDateTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return Date.from(((LocalDate) value).atStartOfDay(ZoneId.systemDefault())
                .toInstant());
    }
}
