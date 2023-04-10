package grok_connect.type.integer;

import grok_connect.type.TypeChecker;

public class DefaultIntegerTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.INTEGER) || (type == java.sql.Types.TINYINT) || (type == java.sql.Types.SMALLINT)
                || typeName.equalsIgnoreCase("int4") || typeName.equalsIgnoreCase("int2")
                || typeName.equalsIgnoreCase("int") || typeName.equalsIgnoreCase("serial2")
                || typeName.equalsIgnoreCase("serial4") || typeName.equalsIgnoreCase("UInt16")
                || typeName.equalsIgnoreCase("UInt8") || (typeName.equalsIgnoreCase("NUMBER")
                && precision < 10 && scale == 0);
    }
}
