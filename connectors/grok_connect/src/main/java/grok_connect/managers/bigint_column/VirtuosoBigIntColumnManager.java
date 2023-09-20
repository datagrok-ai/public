package grok_connect.managers.bigint_column;

import grok_connect.resultset.ColumnMeta;

public class VirtuosoBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return columnMeta.getType() == java.sql.Types.OTHER
                && columnMeta.getPrecision() >= 19 && columnMeta.getScale() == 0;
    }
}
