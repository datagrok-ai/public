package grok_connect.converter.time.impl;

import grok_connect.converter.Converter;
import microsoft.sql.DateTimeOffset;
import java.util.Date;

public class MicrosoftDateTimeOffsetTypeConverter implements Converter<Date> {
    @Override
    public Date convert(Object value) {
        return Date.from(((DateTimeOffset) value).getOffsetDateTime().toInstant());
    }
}
