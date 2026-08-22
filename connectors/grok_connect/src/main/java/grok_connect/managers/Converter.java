package grok_connect.managers;

public interface Converter<T> {
    T convert(Object value);
}
