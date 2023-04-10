package grok_connect.converter.complex.impl;

import com.google.gson.Gson;
import grok_connect.converter.Converter;
import java.util.Map;

public class ComplexTypeConverter implements Converter<Map<String, Object>> {
    private static final Gson gson = new Gson();

    @SuppressWarnings("unchecked")
    @Override
    public Map<String, Object> convert(Object value) {
        if (value instanceof Map) {
            return (Map<String, Object>) value;
        }
        return gson.fromJson(value.toString(), Map.class);
    }
}
