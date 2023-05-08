package grok_connect.column;

import grok_connect.type.DefaultTypeCheckers;
import grok_connect.type.TypeChecker;
import serialization.BigIntColumn;
import serialization.Column;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class BigIntColumnProvider extends AbstractColumnProvider {
    private static final List<TypeChecker> DEFAULT_TYPE_CHECKERS;

    static {
        DEFAULT_TYPE_CHECKERS = new ArrayList<>();
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_BIGINT_TYPECHECKER);
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_BITSTRING_TYPECHECKER);
    }

    public BigIntColumnProvider() {
        super(DEFAULT_TYPE_CHECKERS);
    }

    public BigIntColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new BigIntColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new BigIntColumn(new String[size]);
    }

    @Override
    public boolean isSupported(Object o) {
        return o instanceof Long || o instanceof BigInteger;
    }
}
