package grok_connect.converter.array;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.array.impl.ArrayTypeConverter;
import grok_connect.converter.array.impl.SQLArrayConverter;
import grok_connect.type.TypeChecker;
import oracle.sql.ARRAY;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Array;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class ArrayConverterManager extends AbstractConverterManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ArrayConverterManager.class);
    private static final Converter<String> defaultConverter = Object::toString;
    private final Map<Class<?>, Converter<String>> converterMap;

    public ArrayConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
        converterMap = new HashMap<>();
        SQLArrayConverter sqlArrayConverter = new SQLArrayConverter();
        converterMap.put(Object.class, new ArrayTypeConverter());
        converterMap.put(Array.class, sqlArrayConverter);
        converterMap.put(ARRAY.class, sqlArrayConverter);
    }

    @Override
    public String convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return Arrays.toString(new Object[]{});
        }
        Class<?> aClass = value.getClass();
        Converter<String> converter;
        if (aClass.isArray()) {
            converter = converterMap
                    .get(Object.class);

        } else {
            converter = converterMap
                    .get(aClass);
        }
        if (converter != null) {
            LOGGER.trace("using defined arrayTypeConverter for class {}", aClass);
            return converter.convert(value);
        }
        LOGGER.debug("couldn't find arrayTypeConverter, using default converter for class {}",
                aClass);
        return defaultConverter.convert(value);
    }
}
