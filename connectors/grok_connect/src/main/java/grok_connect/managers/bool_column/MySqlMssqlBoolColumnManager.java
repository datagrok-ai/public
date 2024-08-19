package grok_connect.managers.bool_column;

import grok_connect.resultset.ColumnMeta;

public class MySqlMssqlBoolColumnManager extends DefaultBoolColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        int scale = columnMeta.getScale();
        int precision = columnMeta.getPrecision();
        return type == java.sql.Types.BIT && precision == 1 && scale == 0;
    }
}
