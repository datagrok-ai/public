package grok_connect.managers;

import grok_connect.resultset.ColumnMeta;
import serialization.Column;
import java.sql.ResultSet;
import java.sql.SQLException;

public interface ColumnManager<T> {
    T convert(Object value, ColumnMeta columnMeta);

    boolean isApplicable(ColumnMeta columnMeta);

    Column<?> getColumn(String name, int initColumnSize);

    // resolved once per column in init
    default boolean canReadFast(ColumnMeta columnMeta) {
        return false;
    }

    default void readFast(ResultSet resultSet, int index, Column<?> column, ColumnMeta columnMeta) throws SQLException {
        throw new UnsupportedOperationException();
    }
}
