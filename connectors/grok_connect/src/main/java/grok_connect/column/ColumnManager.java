package grok_connect.column;

import serialization.Column;

public interface ColumnManager<T> {
    T convert(Object value, Object...args);

    boolean isApplicable(int type, String typeName, int precision, int scale);

    default boolean isApplicable(Object o) {
        return false;
    };

    Column getColumn();

    Column getColumnWithInitSize(int size);
}
