package grok_connect.column;

import serialization.Column;

public interface ColumnProvider {
    Column get();

    Column getWithInitSize(int size);

    boolean isSupported(int type, String typeName, int precision, int scale);

    default boolean isSupported(Object o) {
        return false;
    };
}
