package grok_connect.type.time;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DefaultTimeTypeChecker implements TypeChecker {
    private static final Logger LOGGER = LoggerFactory.getLogger(DefaultTimeTypeChecker.class);

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        LOGGER.trace(DEFAULT_LOG_MESSAGE, type, typeName, precision, scale);
        return (type == java.sql.Types.DATE) || (type == java.sql.Types.TIME) || (type == java.sql.Types.TIMESTAMP)
                || type == java.sql.Types.TIMESTAMP_WITH_TIMEZONE
                || type == java.sql.Types.TIME_WITH_TIMEZONE
                || typeName.equalsIgnoreCase("timetz")
                || typeName.equalsIgnoreCase("timestamptz")
                || (typeName.equalsIgnoreCase("TIMESTAMP WITH TIME ZONE"))
                || (typeName.equalsIgnoreCase("datetimeoffset"));
    }
}
