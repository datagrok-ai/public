package grok_connect.managers.integer_column.converters;

import grok_connect.managers.Converter;

public class FloatTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Float) value).intValue();
    }
}
