package grok_connect.converter.float_type.impl;

import grok_connect.converter.Converter;

public class DoubleTypeConverter implements Converter<Float> {
    @Override
    public Float convert(Object value) {
        Double doubleValue = (Double) value;
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
