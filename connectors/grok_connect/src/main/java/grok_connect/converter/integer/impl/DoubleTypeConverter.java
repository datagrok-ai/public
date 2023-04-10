package grok_connect.converter.integer.impl;

import grok_connect.converter.Converter;

public class DoubleTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Double) value).intValue();
    }
}
