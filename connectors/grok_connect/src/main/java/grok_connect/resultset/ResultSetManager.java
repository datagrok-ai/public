package grok_connect.resultset;

import serialization.Column;

public interface ResultSetManager {
    <T> T convert(Object o, int type, String typeName, int precision, int scale, Object...args);

    Column getColumn(int type, String typeName, int precision, int scale);

    Column getColumnWithInitSize(int type, String typeName, int precision, int scale, int size);

    Column getColumn(Object o);
}
