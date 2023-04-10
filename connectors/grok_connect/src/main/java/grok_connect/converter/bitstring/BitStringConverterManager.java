package grok_connect.converter.bitstring;

import grok_connect.converter.AbstractConverterManager;
import grok_connect.converter.Converter;
import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BitStringConverterManager extends AbstractConverterManager<String> {
    private static final Logger LOGGER = LoggerFactory.getLogger(BitStringConverterManager.class);
    private static final Converter<String> defaultConverter = Object::toString;

    public BitStringConverterManager(TypeChecker typeChecker) {
        super(typeChecker);
    }

    @Override
    public String convert(Object value, Object...args) {
        LOGGER.trace("convert method was called");
        LOGGER.trace("using default converter");
        return value == null ? "" : defaultConverter.convert(value);
    }
}
