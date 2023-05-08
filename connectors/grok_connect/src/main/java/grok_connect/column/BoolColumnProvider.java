package grok_connect.column;

import grok_connect.type.DefaultTypeCheckers;
import grok_connect.type.TypeChecker;
import serialization.BoolColumn;
import serialization.Column;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class BoolColumnProvider extends AbstractColumnProvider {
    private static final List<TypeChecker> DEFAULT_TYPE_CHECKERS;

    static {
        DEFAULT_TYPE_CHECKERS = new ArrayList<>();
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_BOOL_TYPECHECKER);
    }

    public BoolColumnProvider() {
        super(DEFAULT_TYPE_CHECKERS);
    }

    public BoolColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new BoolColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new BoolColumn(new Boolean[size]);
    }

    @Override
    public boolean isSupported(Object o) {
        return o instanceof Boolean;
    }
}
