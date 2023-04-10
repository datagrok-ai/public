package grok_connect.converter.bool;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BoolTypeConverterManager extends AbstractConverterManager<Boolean> {
    private static final Logger LOGGER = LoggerFactory.getLogger(BoolTypeConverterManager.class);
    private static final Converter<Boolean> defaultConverter = value -> (Boolean) value;

    public BoolTypeConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
    }

    @Override
    public Boolean convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        if (value == null) {
            LOGGER.trace("value is null");
            return null;
        }
        LOGGER.debug("couldn't find boolTypeConverter, using default converter for class {}",
                value.getClass());
        return defaultConverter.convert(value);
    }
}
