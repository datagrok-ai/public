package grok_connect.type.array;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DefaultArrayTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultArrayTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return type == java.sql.Types.ARRAY || typeName.equalsIgnoreCase("ARRAY");
    }
}
