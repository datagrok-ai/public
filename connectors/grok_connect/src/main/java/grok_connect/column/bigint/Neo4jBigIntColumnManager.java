package grok_connect.column.bigint;

import grok_connect.resultset.ColumnMeta;

import java.sql.Types;

public class Neo4jBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        int type = columnMeta.getType();
        String typeName = columnMeta.getTypeName();
        return typeName.equals("INTEGER") || type == Types.INTEGER;
    }
}
