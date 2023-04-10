package grok_connect.converter;

public interface ConverterManager<T> {
    T convert(Object o, Object...args);

    boolean isSupported(int type, String typeName, int precision, int scale);
}
