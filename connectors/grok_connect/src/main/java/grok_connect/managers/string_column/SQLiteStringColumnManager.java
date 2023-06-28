package grok_connect.managers.string_column;

import grok_connect.resultset.ColumnMeta;
import java.sql.Types;

public class SQLiteStringColumnManager extends DefaultStringColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return columnMeta.getType() == Types.TIMESTAMP || super.isApplicable(columnMeta);
    }
}
