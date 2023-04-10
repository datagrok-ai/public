package grok_connect.type.bigint;

import grok_connect.type.TypeChecker;

import java.sql.Types;

public class Neo4jBigIntTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return typeName.equals("INTEGER") || type == Types.INTEGER;
    }
}
