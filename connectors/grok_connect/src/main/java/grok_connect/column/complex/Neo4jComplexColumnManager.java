package grok_connect.column.complex;

public class Neo4jComplexColumnManager extends DefaultComplexColumnManager {
    @Override
    public boolean isApplicable(int type, String typeName, int precision, int scale) {
        return typeName.equalsIgnoreCase("NODE") || typeName.equalsIgnoreCase("map");
    }
}
