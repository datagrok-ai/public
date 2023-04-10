package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.ComplexTypeColumn;
import java.util.Collection;

public class ComplexTypeColumnProvider extends AbstractColumnProvider {
    public ComplexTypeColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new ComplexTypeColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new ComplexTypeColumn(new Column[size]);
    }
}
