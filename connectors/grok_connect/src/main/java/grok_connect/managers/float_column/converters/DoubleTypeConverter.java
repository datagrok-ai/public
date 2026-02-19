package grok_connect.managers.float_column.converters;

import grok_connect.managers.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DoubleTypeConverter implements Converter<Float> {
    private static final Logger LOGGER = LoggerFactory.getLogger(DoubleTypeConverter.class);

    @Override
    public Float convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        double doubleValue = (Double) value;
        LOGGER.trace("value is {}", doubleValue);
        if (doubleValue > Float.MAX_VALUE) {
            return Float.POSITIVE_INFINITY;
        } else if (doubleValue < - Float.MAX_VALUE) {
            return Float.NEGATIVE_INFINITY;
        } else if (Double.isNaN(doubleValue)) {
            return Float.NaN;
        }
        return (float) doubleValue;
    }
}
