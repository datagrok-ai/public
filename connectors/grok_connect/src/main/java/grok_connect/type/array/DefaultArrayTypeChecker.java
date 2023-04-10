package grok_connect.type.array;

import grok_connect.type.TypeChecker;

public class DefaultArrayTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return type == java.sql.Types.ARRAY || typeName.equalsIgnoreCase("ARRAY");
    }
}
