package grok_connect.converter;

import grok_connect.type.TypeChecker;

public abstract class AbstractConverterManager<T> implements ConverterManager<T> {
    private final TypeChecker typeChecker;

    protected AbstractConverterManager(TypeChecker typeChecker) {
        this.typeChecker = typeChecker;
    }

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return typeChecker.isSupported(type, typeName, precision, scale);
    }
}
