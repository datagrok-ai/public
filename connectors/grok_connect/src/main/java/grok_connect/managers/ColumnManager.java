package grok_connect.managers;

import grok_connect.resultset.ColumnMeta;
import serialization.Column;

public interface ColumnManager<T> {
    T convert(Object value, String columnLabel);

    boolean isApplicable(ColumnMeta columnMeta);

    default boolean isApplicable(Object o) {
        return false;
    }

    Column getColumn();
}
