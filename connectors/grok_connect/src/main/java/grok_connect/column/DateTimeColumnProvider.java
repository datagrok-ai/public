package grok_connect.column;

import grok_connect.type.DefaultTypeCheckers;
import grok_connect.type.TypeChecker;
import microsoft.sql.DateTimeOffset;
import oracle.sql.DATE;
import oracle.sql.TIMESTAMP;
import oracle.sql.TIMESTAMPTZ;
import serialization.Column;
import serialization.DateTimeColumn;
import java.time.temporal.Temporal;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;

public class DateTimeColumnProvider extends AbstractColumnProvider {
    private static final List<TypeChecker> DEFAULT_TYPE_CHECKERS;

    static {
        DEFAULT_TYPE_CHECKERS = new ArrayList<>();
        DEFAULT_TYPE_CHECKERS.add(DefaultTypeCheckers.DEFAULT_DATETIME_TYPECHECKER);
    }

    public DateTimeColumnProvider() {
        super(DEFAULT_TYPE_CHECKERS);
    }

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

    @Override
    public boolean isSupported(Object o) {
        return o instanceof Temporal || o instanceof Date
                || o instanceof DateTimeOffset || o instanceof DATE
                || o instanceof TIMESTAMP || o instanceof TIMESTAMPTZ;
    }
}
