package grok_connect.managers.bigint_column;

import grok_connect.resultset.ColumnMeta;
import java.sql.Types;

public class ClickHouseBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        return type == Types.BIGINT
                || typeName.equalsIgnoreCase("UInt32")
                || typeName.equalsIgnoreCase("UInt64")
                || typeName.equalsIgnoreCase("UInt128")
                || typeName.equalsIgnoreCase("UInt256")
                || typeName.equalsIgnoreCase("Int64")
                || typeName.equalsIgnoreCase("Int128")
                || typeName.equalsIgnoreCase("Int256");
    }
}
