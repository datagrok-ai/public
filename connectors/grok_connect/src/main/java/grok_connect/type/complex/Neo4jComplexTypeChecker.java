package grok_connect.type.complex;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Neo4jComplexTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(Neo4jComplexTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return typeName.equalsIgnoreCase("NODE");
    }
}
