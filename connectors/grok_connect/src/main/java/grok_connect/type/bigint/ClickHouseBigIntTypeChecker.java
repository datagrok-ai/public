package grok_connect.type.bigint;

import grok_connect.type.TypeChecker;

import java.sql.Types;

public class ClickHouseBigIntTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
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
