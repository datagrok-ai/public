package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.FloatColumn;

import java.util.Collection;

public class FloatColumnProvider extends AbstractColumnProvider {
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
}
