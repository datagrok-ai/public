package grok_connect.column.integer;

public class Neo4jIntColumnManager extends DefaultIntColumnManager {
    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return false;
    }
}
