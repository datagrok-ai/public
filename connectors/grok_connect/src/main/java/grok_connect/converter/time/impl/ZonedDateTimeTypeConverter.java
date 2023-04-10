package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import java.time.ZonedDateTime;
import java.util.Date;

public class ZonedDateTimeTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return Date.from(((ZonedDateTime)value).toInstant());
    }
}
