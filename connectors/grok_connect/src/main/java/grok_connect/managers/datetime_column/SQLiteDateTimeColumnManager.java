package grok_connect.managers.datetime_column;

import grok_connect.resultset.ColumnMeta;

public class SQLiteDateTimeColumnManager extends DefaultDateTimeColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return false;
    }
}
