package grok_connect.managers.integer_column;

import grok_connect.resultset.ColumnMeta;

public class Neo4jIntColumnManager extends DefaultIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        return false;
    }
}
