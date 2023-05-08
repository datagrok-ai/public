package grok_connect.column;

import grok_connect.type.DefaultTypeCheckers;
import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.IntColumn;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class IntColumnProvider extends AbstractColumnProvider {
    private static final List<TypeChecker> DEFAULT_TYPE_CHECKERS;

    static {
        DEFAULT_TYPE_CHECKERS = new ArrayList<>();
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_INT_TYPECHECKER);
    }

    public IntColumnProvider() {
        super(DEFAULT_TYPE_CHECKERS);
    }

    public IntColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new IntColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new IntColumn(new Integer[size]);
    }

    @Override
    public boolean isSupported(Object o) {
        return o instanceof Byte || o instanceof Short || o instanceof Integer;
    }
}
