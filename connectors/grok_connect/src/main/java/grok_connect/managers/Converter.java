package grok_connect.managers;

public interface Converter<T> {
    String DEFAULT_LOG_MESSAGE = "convert method was called for object with class {}";

    T convert(Object value);
}
