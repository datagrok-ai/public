package grok_connect.type.integer;

import grok_connect.type.TypeChecker;

public class OracleSnowflakeIntegerTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return typeName.equalsIgnoreCase("number") && precision < 10 && scale == 0;
    }
}
