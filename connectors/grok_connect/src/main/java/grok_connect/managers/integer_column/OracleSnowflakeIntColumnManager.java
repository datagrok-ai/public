package grok_connect.managers.integer_column;

import grok_connect.resultset.ColumnMeta;

public class OracleSnowflakeIntColumnManager extends DefaultIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        String typeName = columnMeta.getTypeName();
        int precision = columnMeta.getPrecision();
        int scale = columnMeta.getScale();
        return typeName.equalsIgnoreCase("number") && precision < 10 && scale == 0;
    }
}
