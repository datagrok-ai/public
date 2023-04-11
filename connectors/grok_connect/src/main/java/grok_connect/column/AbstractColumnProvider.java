package grok_connect.column;

import grok_connect.type.TypeChecker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.util.Collection;

public abstract class AbstractColumnProvider implements ColumnProvider {
    protected Logger logger = LoggerFactory.getLogger(this.getClass());
    private final Collection<TypeChecker> typeCheckers;

    protected AbstractColumnProvider(Collection<TypeChecker> typeCheckers) {
        this.typeCheckers = typeCheckers;
    }

    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        logger.trace("method isSupported was called with arguments: type {}, typeName {}, "
                + "precision {}, scale {}", type, typeName, precision, scale);
        return typeCheckers.stream()
                .anyMatch(checker -> checker.isSupported(type, typeName, precision, scale));
    }
}
