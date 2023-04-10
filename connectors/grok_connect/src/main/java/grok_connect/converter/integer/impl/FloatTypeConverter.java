package grok_connect.converter.integer.impl;

import grok_connect.converter.Converter;

public class FloatTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Float) value).intValue();
    }
}
