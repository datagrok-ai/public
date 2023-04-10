package grok_connect.type.bitstring;

import grok_connect.type.TypeChecker;

public class DefaultBitStringTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.BIT && precision > 1 || scale > 1 && type == java.sql.Types.BIT)
                || typeName.equalsIgnoreCase("varbit");
    }
}
