package grok_connect.type.bigint;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.sql.Types;

public class ClickHouseBigIntTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(ClickHouseBigIntTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return type == Types.BIGINT
                || typeName.equalsIgnoreCase("UInt32")
                || typeName.equalsIgnoreCase("UInt64")
                || typeName.equalsIgnoreCase("UInt128")
                || typeName.equalsIgnoreCase("UInt256")
                || typeName.equalsIgnoreCase("Int64")
                || typeName.equalsIgnoreCase("Int128")
                || typeName.equalsIgnoreCase("Int256");
    }
}
