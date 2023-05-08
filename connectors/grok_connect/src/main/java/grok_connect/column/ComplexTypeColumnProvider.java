package grok_connect.column;

import grok_connect.type.DefaultTypeCheckers;
import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.ComplexTypeColumn;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

public class ComplexTypeColumnProvider extends AbstractColumnProvider {
    private static final List<TypeChecker> DEFAULT_TYPE_CHECKERS;

    static {
        DEFAULT_TYPE_CHECKERS = new ArrayList<>();
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_COMPLEX_TYPECHECKER);
    }

    public ComplexTypeColumnProvider() {
        super(DEFAULT_TYPE_CHECKERS);
    }

    public ComplexTypeColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new ComplexTypeColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new ComplexTypeColumn(new Column[size]);
    }

    @Override
    public boolean isSupported(Object o) {
        return o instanceof Map;
    }
}
