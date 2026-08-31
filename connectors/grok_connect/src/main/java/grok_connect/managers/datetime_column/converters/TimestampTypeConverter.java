package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import java.sql.Timestamp;
import java.util.Date;

public class TimestampTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return (Timestamp) value;
    }
}
