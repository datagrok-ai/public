package grok_connect.converter.float_type.impl;

import grok_connect.converter.Converter;
import java.math.BigDecimal;

public class BigDecimalTypeConverter implements Converter<Float> {
    @Override
    public Float convert(Object value) {
        return ((BigDecimal)value).floatValue();
    }
}
