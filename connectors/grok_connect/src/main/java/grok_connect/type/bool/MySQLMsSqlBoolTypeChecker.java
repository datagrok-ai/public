package grok_connect.type.bool;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Types;

public class MySQLMsSqlBoolTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(MySQLMsSqlBoolTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return type == Types.BIT && precision == 1 && scale == 0;
    }
}
