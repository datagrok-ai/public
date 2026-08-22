package grok_connect.managers.integer_column.converters;

import grok_connect.managers.Converter;

public class LongTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Long) value).intValue();
    }
}
