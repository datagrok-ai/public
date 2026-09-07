package grok_connect.managers.bool_column;

import grok_connect.managers.ColumnManager;
import grok_connect.managers.Converter;
import grok_connect.resultset.ColumnMeta;
import serialization.BoolColumn;
import serialization.Column;
import java.sql.ResultSet;
import java.sql.SQLException;

public class DefaultBoolColumnManager implements ColumnManager<Boolean> {
    private static final Converter<Boolean> DEFAULT_CONVERTER = value -> (Boolean) value;

    @Override
    public Boolean convert(Object value, ColumnMeta columnMeta) {
        if (value == null)
            return null;
        return DEFAULT_CONVERTER.convert(value);
    }

    @Override
    public boolean canReadFast(ColumnMeta columnMeta) {
        return isApplicable(columnMeta);
    }

    // SQL NULL reads as false, exactly what BoolColumn.add(null) stores.
    @Override
    public void readFast(ResultSet resultSet, int index, Column<?> column, ColumnMeta columnMeta) throws SQLException {
        ((BoolColumn) column).add(resultSet.getBoolean(index));
    }

    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        int precision = columnMeta.getPrecision();
        int scale = columnMeta.getScale();
        return (type == java.sql.Types.BOOLEAN) ||
                typeName.equalsIgnoreCase("bool") || (type == java.sql.Types.BIT
                && precision == 1 && scale == 0);
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new BoolColumn(name, initColumnSize);
    }
}
