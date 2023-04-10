package grok_connect.type.bool;

import grok_connect.type.TypeChecker;

public class DefaultBoolTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.BOOLEAN) ||
                typeName.equalsIgnoreCase("bool") || (type == java.sql.Types.BIT
                && scale == 1 && precision == 0);
    }
}
