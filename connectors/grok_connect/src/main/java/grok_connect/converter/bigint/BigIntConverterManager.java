package grok_connect.converter.bigint;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BigIntConverterManager extends AbstractConverterManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(BigIntConverterManager.class);
    private static final Converter<String> defaultConverter = Object::toString;

    public BigIntConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
    }

    @Override
    public String convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        LOGGER.trace("using default converter");
        return value == null ? "" : defaultConverter.convert(value);
    }
}
