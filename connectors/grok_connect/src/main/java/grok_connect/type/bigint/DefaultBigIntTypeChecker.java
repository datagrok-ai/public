package grok_connect.type.bigint;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DefaultBigIntTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultBigIntTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return (type == java.sql.Types.BIGINT) || typeName.equalsIgnoreCase("int8") ||
                typeName.equalsIgnoreCase("serial8");
    }
}
