package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.IntColumn;
import java.util.Collection;

public class IntColumnProvider extends AbstractColumnProvider {
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
}
