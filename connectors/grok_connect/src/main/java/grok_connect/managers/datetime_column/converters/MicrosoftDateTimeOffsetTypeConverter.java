package grok_connect.managers.datetime_column.converters;

import grok_connect.managers.Converter;
import microsoft.sql.DateTimeOffset;
import java.sql.Timestamp;
import java.util.Date;

public class MicrosoftDateTimeOffsetTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        // Timestamp keeps the sub-ms fraction (Date.from would truncate to ms)
        return Timestamp.from(((DateTimeOffset) value).getOffsetDateTime().toInstant());
    }
}
