package grok_connect.column.bool;

public class MySqlMssqlBoolColumnManager extends DefaultBoolColumnManager {
    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return type == java.sql.Types.BIT && precision == 1 && scale == 0;
    }
}
