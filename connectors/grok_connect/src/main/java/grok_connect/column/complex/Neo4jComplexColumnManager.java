package grok_connect.column.complex;

import grok_connect.resultset.ColumnMeta;

public class Neo4jComplexColumnManager extends DefaultComplexColumnManager {
    @Override
    public boolean isApplicable(ColumnMeta columnMeta) {
        String typeName = columnMeta.getTypeName();
        return typeName.equalsIgnoreCase("NODE") || typeName.equalsIgnoreCase("map");
    }
}
