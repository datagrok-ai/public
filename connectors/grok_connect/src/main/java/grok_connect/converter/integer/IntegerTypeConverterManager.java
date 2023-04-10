package grok_connect.converter.integer;

import com.clickhouse.data.value.UnsignedByte;
import com.clickhouse.data.value.UnsignedShort;
import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.integer.impl.*;
import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;

public class IntegerTypeConverterManager extends AbstractConverterManager<Integer> {
    private static final Logger LOGGER = LoggerFactory.getLogger(IntegerTypeConverterManager.class);
    private static final Converter<Integer> defaultConverter = value -> (Integer) value;
    private final Map<Class<?>, Converter<Integer>> converterMap;

    public IntegerTypeConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
        converterMap = new HashMap<>();
        converterMap.put(BigDecimal.class, new BigDecimalTypeConverter());
        converterMap.put(Byte.class, new ByteTypeConverter());
        converterMap.put(UnsignedByte.class, new ClickHouseUnsignedByteTypeConverter());
        converterMap.put(UnsignedShort.class, new ClickHouseUnsignedShortTypeConverter());
        converterMap.put(Double.class, new DoubleTypeConverter());
        converterMap.put(Float.class, new FloatTypeConverter());
        converterMap.put(Long.class, new LongTypeConverter());
        converterMap.put(Short.class, new ShortTypeConverter());
    }

    @Override
    public Integer convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<Integer> converter = converterMap
                .get(aClass);
        if (converter != null) {
            LOGGER.trace("using defined integerTypeConverter for class {}", aClass);
            return converter.convert(value);
        }
        LOGGER.debug("couldn't find integerTypeConverter, using default converter for class {}",
                aClass);
        return defaultConverter.convert(value);
    }
}
