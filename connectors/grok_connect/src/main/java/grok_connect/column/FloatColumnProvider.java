package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.FloatColumn;

import java.math.BigDecimal;
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

    @Override
    public boolean isSupported(Object o) {
        return o instanceof Float || o instanceof Double || o instanceof BigDecimal;
    }
}
