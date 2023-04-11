package grok_connect.converter.float_type.impl;

import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DoubleTypeConverter implements Converter<Float> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DoubleTypeConverter.class);

    @Override
    public Float convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        Double doubleValue = (Double) value;
        LOGGER.trace("value is {}", doubleValue);
        if (doubleValue == Double.POSITIVE_INFINITY || doubleValue > Float.MAX_VALUE) {
            return Float.POSITIVE_INFINITY;
        } else if (doubleValue == Double.NEGATIVE_INFINITY || doubleValue < - Float.MAX_VALUE) {
            return Float.NEGATIVE_INFINITY;
        } else if (Double.isNaN(doubleValue)) {
            return Float.NaN;
        }
        return doubleValue.floatValue();
    }
}
