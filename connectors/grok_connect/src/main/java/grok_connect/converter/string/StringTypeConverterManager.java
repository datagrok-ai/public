package grok_connect.converter.string;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.converter.string.impl.ClobTypeConverter;
import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Clob;
import java.util.HashMap;
import java.util.Map;
import com.ibm.db2.jcc.am.c9;

public class StringTypeConverterManager extends AbstractConverterManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(StringTypeConverterManager.class);
    private static final Converter<String> defaultConverter = Object::toString;
    private final Map<Class<?>, Converter<String>> converterMap;

    public StringTypeConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
        converterMap = new HashMap<>();
        Converter<String> clobConverter = new ClobTypeConverter();
        converterMap.put(Clob.class, clobConverter);
        converterMap.put(c9.class, clobConverter);
    }

    @Override
    public String convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        Class<?> aClass = value.getClass();
        Converter<String> converter = converterMap
                .get(aClass);
        if (converter != null) {
            LOGGER.trace("using defined stringTypeConverter for class {}", aClass);
            return converter.convert(value);
        }
        LOGGER.debug("couldn't find stringTypeConverter, using default converter for class {}",
                aClass);
        return defaultConverter.convert(value);
    }
}
