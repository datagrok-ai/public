package grok_connect.managers.datetime_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.managers.datetime_column.converters.*;
import grok_connect.resultset.ColumnMeta;
import microsoft.sql.DateTimeOffset;
import oracle.sql.TIMESTAMPTZ;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import serialization.Column;
import serialization.DateTimeColumn;
import java.sql.Timestamp;
import java.time.*;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

public class DefaultDateTimeColumnManager implements ColumnManager<Double> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultDateTimeColumnManager.class);
    private static final Converter<Date> DEFAULT_CONVERTER = value -> (Date) value;
    private static final double DEFAULT_MULTIPLIER = 1000.0;
    private final Map<Class<?>, Converter<Date>> converterMap;

    {
        converterMap = new HashMap<>();
        converterMap.put(LocalDateTime.class, new LocalDateTimeTypeConverter());
        converterMap.put(LocalDate.class, new LocalDateTypeConverter());
        converterMap.put(DateTimeOffset.class, new MicrosoftDateTimeOffsetTypeConverter());
        converterMap.put(OffsetDateTime.class, new OffsetDateTimeTypeConverter());
        converterMap.put(TIMESTAMPTZ.class, new OracleTimestampTZTypeConverter());
        converterMap.put(Timestamp.class, new TimestampTypeConverter());
        converterMap.put(ZonedDateTime.class, new ZonedDateTimeTypeConverter());
        converterMap.put(Instant.class, new InstantTypeConverter());
        converterMap.put(LocalTime.class, new LocalTimeConverter());
    }

    @Override
    public Double convert(Object value, ColumnMeta columnMeta) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<Date> converter = converterMap
                .getOrDefault(aClass, DEFAULT_CONVERTER);
        Date date =  converter.convert(value);
        return date.getTime() * DEFAULT_MULTIPLIER;
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        return (type == java.sql.Types.DATE) || (type == java.sql.Types.TIME) ||
                (type == java.sql.Types.TIMESTAMP)
                || type == java.sql.Types.TIMESTAMP_WITH_TIMEZONE
                || type == java.sql.Types.TIME_WITH_TIMEZONE
                || typeName.equalsIgnoreCase("timetz")
                || typeName.equalsIgnoreCase("timestamptz")
                || (typeName.equalsIgnoreCase("TIMESTAMP WITH TIME ZONE"))
                || (typeName.equalsIgnoreCase("datetimeoffset"));
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new DateTimeColumn(name, initColumnSize);
    }
}
