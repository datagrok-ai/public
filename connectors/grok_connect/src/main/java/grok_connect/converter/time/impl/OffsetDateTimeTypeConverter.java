package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import java.time.OffsetDateTime;
import java.util.Date;

public class OffsetDateTimeTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return Date.from(((OffsetDateTime) value)
                .toInstant());
    }
}
