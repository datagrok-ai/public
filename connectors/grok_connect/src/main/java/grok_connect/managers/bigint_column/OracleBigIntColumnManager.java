package grok_connect.managers.bigint_column;

import grok_connect.resultset.ColumnMeta;
import serialization.Column;
import serialization.IntColumn;

public class OracleBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        String typeName = columnMeta.getTypeName();
        int precision = columnMeta.getPrecision();
        int scale = columnMeta.getScale();
        return typeName.equalsIgnoreCase("number") && precision > 10 && scale == 0;
    }

    @Override
    public Column<?> getColumn(String name, int initColumnSize) {
        return new IntColumn(name, initColumnSize);
    }
}
