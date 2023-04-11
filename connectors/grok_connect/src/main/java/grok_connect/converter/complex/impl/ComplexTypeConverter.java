package grok_connect.converter.complex.impl;

import com.google.gson.Gson;
import grok_connect.converter.Converter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.Map;

public class ComplexTypeConverter implements Converter<Map<String, Object>> {
    private static final Logger LOGGER = LoggerFactory.getLogger(ComplexTypeConverter.class);
    private static final Gson gson = new Gson();

    @SuppressWarnings("unchecked")
    @Override
    public Map<String, Object> convert(Object value) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, value.getClass());
        if (value instanceof Map) {
            LOGGER.trace("value is map, converting it to Map<String, Object>");
            return (Map<String, Object>) value;
        }
        LOGGER.trace("value is not map, trying to convert string representation to map");
        return gson.fromJson(value.toString(), Map.class);
    }
}
