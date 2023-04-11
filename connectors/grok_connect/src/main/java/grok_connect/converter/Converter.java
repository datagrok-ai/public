package grok_connect.converter;

import java.util.Date;

public interface Converter<T> {
    String DEFAULT_LOG_MESSAGE = "convert method was called with for object with class {}";

    T convert(Object value);
}
