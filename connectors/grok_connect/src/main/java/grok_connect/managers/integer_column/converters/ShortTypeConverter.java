package grok_connect.managers.integer_column.converters;

import grok_connect.managers.Converter;

public class ShortTypeConverter implements Converter<Integer> {
    @Override
    public Integer convert(Object value) {
        return ((Short)value).intValue();
    }
}
