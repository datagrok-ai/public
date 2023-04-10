package grok_connect.type.float_type;

import grok_connect.type.TypeChecker;

import java.sql.Types;

public class DefaultFloatTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return type == Types.FLOAT || type == java.sql.Types.DOUBLE || type == java.sql.Types.REAL ||
                type == Types.DECIMAL ||
                typeName.equalsIgnoreCase("float8") ||
                typeName.equalsIgnoreCase("float4") ||
                typeName.equalsIgnoreCase("money") ||
                typeName.equalsIgnoreCase("binary_float") ||
                typeName.equalsIgnoreCase("binary_double") ||
                typeName.equalsIgnoreCase("numeric") ||
                typeName.equalsIgnoreCase("DECFLOAT") ||
                (typeName.equalsIgnoreCase("number") && scale > 0);
    }
}
