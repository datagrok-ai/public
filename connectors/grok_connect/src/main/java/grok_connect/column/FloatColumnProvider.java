package grok_connect.column;

import grok_connect.type.DefaultTypeCheckers;
import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.FloatColumn;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class FloatColumnProvider extends AbstractColumnProvider {
    private static final List<TypeChecker> DEFAULT_TYPE_CHECKERS;

    static {
        DEFAULT_TYPE_CHECKERS = new ArrayList<>();
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_FLOAT_TYPECHECKER);
    }

    public FloatColumnProvider() {
        super(DEFAULT_TYPE_CHECKERS);
    }

    public FloatColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new FloatColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new FloatColumn(new Float[size]);
    }

    @Override
    public boolean isSupported(Object o) {
        return o instanceof Float || o instanceof Double || o instanceof BigDecimal;
    }
}
