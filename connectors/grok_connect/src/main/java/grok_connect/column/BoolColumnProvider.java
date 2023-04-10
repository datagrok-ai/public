package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.BoolColumn;
import serialization.Column;
import java.util.Collection;

public class BoolColumnProvider extends AbstractColumnProvider {
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
}
