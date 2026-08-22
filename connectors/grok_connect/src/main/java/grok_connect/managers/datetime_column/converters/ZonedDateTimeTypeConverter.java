package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import java.sql.Timestamp;
import java.time.ZonedDateTime;
import java.util.Date;

public class ZonedDateTimeTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        // Timestamp keeps the sub-ms fraction (Date.from would truncate to ms)
        return Timestamp.from(((ZonedDateTime)value).toInstant());
    }
}
