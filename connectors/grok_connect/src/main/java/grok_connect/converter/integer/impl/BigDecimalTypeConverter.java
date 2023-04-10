package grok_connect.converter.integer.impl;

import grok_connect.converter.Converter;

import java.math.BigDecimal;

public class BigDecimalTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((BigDecimal) value).intValue();
    }
}
