package grok_connect.type.time;

import grok_connect.type.TypeChecker;

public class DefaultTimeTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.DATE) || (type == java.sql.Types.TIME) || (type == java.sql.Types.TIMESTAMP)
                || type == java.sql.Types.TIMESTAMP_WITH_TIMEZONE
                || type == java.sql.Types.TIME_WITH_TIMEZONE
                || typeName.equalsIgnoreCase("timetz")
                || typeName.equalsIgnoreCase("timestamptz")
                || (typeName.equalsIgnoreCase("TIMESTAMP WITH TIME ZONE"))
                || (typeName.equalsIgnoreCase("datetimeoffset"));
    }
}
