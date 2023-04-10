package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.BigIntColumn;
import serialization.Column;

import java.util.Collection;

public class BigIntColumnProvider extends AbstractColumnProvider {
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
}
