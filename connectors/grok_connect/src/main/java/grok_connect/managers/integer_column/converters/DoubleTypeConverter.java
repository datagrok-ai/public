package grok_connect.managers.integer_column.converters;

import grok_connect.managers.Converter;

public class DoubleTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Double) value).intValue();
    }
}
