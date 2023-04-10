package grok_connect.column;

import grok_connect.type.TypeChecker;
import serialization.Column;
import serialization.DateTimeColumn;
import java.util.Collection;

public class DateTimeColumnProvider extends AbstractColumnProvider {
    public DateTimeColumnProvider(Collection<TypeChecker> typeCheckers) {
        super(typeCheckers);
    }

    @Override
    public Column get() {
        return new DateTimeColumn();
    }

    @Override
    public Column getWithInitSize(int size) {
        return new DateTimeColumn(new Double[size]);
    }
}
