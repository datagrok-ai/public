package grok_connect.managers.float_column.converters;

import grok_connect.managers.Converter;

public class DoubleTypeConverter implements Converter<Float> {
    @Override
    public Float convert(Object value) {
        return narrow((Double) value);
    }

    public static float narrow(double doubleValue) {
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
