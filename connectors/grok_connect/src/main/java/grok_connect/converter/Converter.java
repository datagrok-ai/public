package grok_connect.converter;

import java.util.Date;

public interface Converter<T> {
    T convert(Object value);
}
