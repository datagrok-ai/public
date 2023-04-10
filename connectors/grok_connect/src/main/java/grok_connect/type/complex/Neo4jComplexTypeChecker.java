package grok_connect.type.complex;

import grok_connect.type.TypeChecker;

public class Neo4jComplexTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return typeName.equalsIgnoreCase("NODE");
    }
}
