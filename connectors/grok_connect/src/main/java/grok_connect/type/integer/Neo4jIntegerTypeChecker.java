package grok_connect.type.integer;

import grok_connect.type.TypeChecker;

public class Neo4jIntegerTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return false;
    }
}
