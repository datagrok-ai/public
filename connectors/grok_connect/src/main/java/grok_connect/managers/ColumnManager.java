package grok_connect.managers;

import grok_connect.resultset.ColumnMeta;
import serialization.Column;

public interface ColumnManager<T> {
    T convert(Object value, ColumnMeta columnMeta);

    boolean isApplicable(ColumnMeta columnMeta);

    Column<?> getColumn(String name, int initColumnSize);
}
