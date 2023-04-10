package grok_connect.column;

import grok_connect.type.TypeChecker;

import java.util.Collection;

public abstract class AbstractColumnProvider implements ColumnProvider {
    private final Collection<TypeChecker> typeCheckers;

    protected AbstractColumnProvider(Collection<TypeChecker> typeCheckers) {
        this.typeCheckers = typeCheckers;
    }

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return typeCheckers.stream()
                .anyMatch(checker -> checker.isSupported(type, typeName, precision, scale));
    }
}
