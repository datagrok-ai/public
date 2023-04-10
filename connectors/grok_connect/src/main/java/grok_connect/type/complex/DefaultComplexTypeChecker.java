package grok_connect.type.complex;

import grok_connect.type.TypeChecker;

public class DefaultComplexTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return false;
    }
}
