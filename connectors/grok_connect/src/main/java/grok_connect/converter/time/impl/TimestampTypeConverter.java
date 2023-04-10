package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import java.sql.Timestamp;
import java.util.Date;

public class TimestampTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return (Timestamp) value;
    }
}
