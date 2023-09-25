package grok_connect.managers.bigint_column;

import grok_connect.resultset.ColumnMeta;

public class SnowflakeBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        String typeName = columnMeta.getTypeName();
        int precision = columnMeta.getPrecision();
        int scale = columnMeta.getScale();
        return typeName.equalsIgnoreCase("number") && precision > 10 && scale == 0;
    }
}
