package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.StringColumn;

import java.util.Collection;

public class StringColumnProvider extends AbstractColumnProvider {
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
