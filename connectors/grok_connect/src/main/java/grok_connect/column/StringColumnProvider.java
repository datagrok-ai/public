package grok_connect.column;

import grok_connect.type.DefaultTypeCheckers;
import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.StringColumn;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class StringColumnProvider extends AbstractColumnProvider {
    private static final List<TypeChecker> DEFAULT_TYPE_CHECKERS;

    static {
        DEFAULT_TYPE_CHECKERS = new ArrayList<>();
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_ARRAY_TYPECHECKER);
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_STRING_TYPECHECKER);
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_XML_TYPECHECKER);
    }

    public StringColumnProvider() {
        super(DEFAULT_TYPE_CHECKERS);
    }

    public StringColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new StringColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new StringColumn(new String[size]);
    }

    @Override
    public boolean isSupported(Object o) {
        return o instanceof String;
    }
}
