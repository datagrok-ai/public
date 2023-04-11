package grok_connect.type.string;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DefaultStringTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultStringTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return ((type == java.sql.Types.VARCHAR)|| (type == java.sql.Types.CHAR) ||
                (type == java.sql.Types.LONGVARCHAR) || (type == java.sql.Types.CLOB)
                || (type == java.sql.Types.NCLOB) ||
                typeName.equalsIgnoreCase("varchar") ||
                typeName.equalsIgnoreCase("nvarchar") ||
                typeName.equalsIgnoreCase("nchar") ||
                typeName.equalsIgnoreCase("ntext")) &&
                !typeName.equalsIgnoreCase("uuid") &&
                !typeName.equalsIgnoreCase("set");
    }
}
