package grok_connect.column.bigint;

import java.sql.Types;

public class Neo4jBigIntColumnManager extends DefaultBigIntColumnManager {
    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return typeName.equals("INTEGER") || type == Types.INTEGER;
    }
}
