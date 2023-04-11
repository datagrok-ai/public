package grok_connect.converter;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
