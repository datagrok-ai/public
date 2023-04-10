package grok_connect.converter.float_type;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.float_type.impl.BigDecimalTypeConverter;
import grok_connect.converter.float_type.impl.DoubleTypeConverter;
import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;

public class FloatTypeConverterManager extends AbstractConverterManager<Float> {
    private static final Logger LOGGER = LoggerFactory.getLogger(FloatTypeConverterManager.class);
    private static final Converter<Float> defaultConverter = value -> (Float) value;
    private final Map<Class<?>, Converter<Float>> converterMap;

    public FloatTypeConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
        converterMap = new HashMap<>();
        converterMap.put(Double.class, new DoubleTypeConverter());
        converterMap.put(BigDecimal.class, new BigDecimalTypeConverter());
    }

    @Override
    public Float convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<Float> converter = converterMap
                .get(aClass);
        if (converter != null) {
            LOGGER.trace("using defined floatTypeConverter for class {}", aClass);
            return converter.convert(value);
        }
        LOGGER.debug("couldn't find floatTypeConverter, using default converter for class {}",
                aClass);
        return defaultConverter.convert(value);
    }
}
