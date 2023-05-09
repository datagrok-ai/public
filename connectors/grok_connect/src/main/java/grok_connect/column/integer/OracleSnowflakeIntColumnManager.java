package grok_connect.column.integer;

public class OracleSnowflakeIntColumnManager extends DefaultIntColumnManager {
    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return typeName.equalsIgnoreCase("number") && precision < 10 && scale == 0;
    }
}
