package grok_connect.type.bigint;

import grok_connect.type.TypeChecker;

public class DefaultBigIntTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.BIGINT) || typeName.equalsIgnoreCase("int8") ||
                typeName.equalsIgnoreCase("serial8");
    }
}
