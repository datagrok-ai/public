package grok_connect.column.bigint;

import java.sql.Types;

public class ClickHouseBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
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
