package grok_connect.converter.time;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.time.impl.LocalDateTimeTypeConverter;
import grok_connect.converter.time.impl.LocalDateTypeConverter;
import grok_connect.converter.time.impl.MicrosoftDateTimeOffsetTypeConverter;
import grok_connect.converter.time.impl.OffsetDateTimeTypeConverter;
import grok_connect.converter.time.impl.OracleTimestampTZTypeConverter;
import grok_connect.converter.time.impl.TimestampTypeConverter;
import grok_connect.converter.time.impl.ZonedDateTimeTypeConverter;
import grok_connect.type.TypeChecker;
import microsoft.sql.DateTimeOffset;
import oracle.sql.TIMESTAMPTZ;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Timestamp;
import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.OffsetDateTime;
import java.time.ZonedDateTime;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

public class TimeTypeConverterManager extends AbstractConverterManager<Double> {
    private static final Logger LOGGER = LoggerFactory.getLogger(TimeTypeConverterManager.class);
    private static final Converter<Date> defaultConverter = value -> (Date) value;
    private static final double DEFAULT_MULTIPLIER = 1000.0;
    private final Map<Class<?>, Converter<Date>> converterMap;

    public TimeTypeConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
        converterMap = new HashMap<>();
        converterMap.put(LocalDateTime.class, new LocalDateTimeTypeConverter());
        converterMap.put(LocalDate.class, new LocalDateTypeConverter());
        converterMap.put(DateTimeOffset.class, new MicrosoftDateTimeOffsetTypeConverter());
        converterMap.put(OffsetDateTime.class, new OffsetDateTimeTypeConverter());
        converterMap.put(TIMESTAMPTZ.class, new OracleTimestampTZTypeConverter());
        converterMap.put(Timestamp.class, new TimestampTypeConverter());
        converterMap.put(ZonedDateTime.class, new ZonedDateTimeTypeConverter());
    }

    @Override
    public Double convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<Date> timeTypeConverter = converterMap
                .get(aClass);
        Date date;
        if (timeTypeConverter != null) {
            LOGGER.trace("using defined timeTypeConverter for class {}", aClass);
            date =  timeTypeConverter.convert(value);
        } else {
            LOGGER.debug("couldn't find timeTypeConverter, using default converter for class {}",
                    aClass);
            date =  defaultConverter.convert(value);
        }
        return date.getTime() * DEFAULT_MULTIPLIER;
    }
}
