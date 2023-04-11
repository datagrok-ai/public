package grok_connect.type.bitstring;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DefaultBitStringTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultBitStringTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return (type == java.sql.Types.BIT && precision > 1 || scale > 1 && type == java.sql.Types.BIT)
                || typeName.equalsIgnoreCase("varbit");
    }
}
