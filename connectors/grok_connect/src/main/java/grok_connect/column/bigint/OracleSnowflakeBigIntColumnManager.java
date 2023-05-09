package grok_connect.column.bigint;

public class OracleSnowflakeBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return typeName.equalsIgnoreCase("number") && precision > 10 && scale == 0;
    }
}
