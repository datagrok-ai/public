package grok_connect.type.float_type;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Types;

public class DefaultFloatTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultFloatTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
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
