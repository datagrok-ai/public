package grok_connect.managers.float_column.converters;

import grok_connect.managers.Converter;

import java.math.BigDecimal;

public class BigDecimalTypeConverter implements Converter<Float> {
    @Override
    public Float convert(Object value) {
        return ((BigDecimal)value).floatValue();
    }
}
