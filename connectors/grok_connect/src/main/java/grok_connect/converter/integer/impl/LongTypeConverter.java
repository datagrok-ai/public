package grok_connect.converter.integer.impl;

import grok_connect.converter.Converter;

public class LongTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Long) value).intValue();
    }
}
