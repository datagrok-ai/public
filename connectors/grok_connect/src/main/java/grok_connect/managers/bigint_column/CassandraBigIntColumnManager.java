package grok_connect.managers.bigint_column;

import grok_connect.resultset.ColumnMeta;

public class CassandraBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return (columnMeta.getType() == java.sql.Types.BIGINT)
                || columnMeta.getTypeName().equalsIgnoreCase("varint");
    }
}
