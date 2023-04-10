package grok_connect.converter.integer.impl;

import grok_connect.converter.Converter;

public class ShortTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Short)value).intValue();
    }
}
