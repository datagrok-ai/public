package grok_connect.managers.string_column.converters;

import grok_connect.managers.Converter;
import java.lang.reflect.Array;

public class ArrayTypeConverter implements Converter<String> {
    @Override
    public String convert(Object value) {
        return getStringArrayRepresentation(value);
    }

    private String getStringArrayRepresentation(Object array) {
        int length = Array.getLength(array);
        StringBuilder builder = new StringBuilder("[");
        for (int i = 0; i < length; i++) {
            Object o = Array.get(array, i);
            if (o.getClass().isArray()) {
                builder.append(getStringArrayRepresentation(o));
            } else {
                builder.append(o);
            }
            if (i != length - 1) {
                builder.append(", ");
            }
        }
        builder.append("]");
        return builder.toString();
    }
}
