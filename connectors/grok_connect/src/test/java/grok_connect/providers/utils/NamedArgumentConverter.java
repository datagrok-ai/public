package grok_connect.providers.utils;

import org.junit.jupiter.api.Named;
import org.junit.jupiter.api.extension.ParameterContext;
import org.junit.jupiter.params.converter.ArgumentConversionException;
import org.junit.jupiter.params.converter.ArgumentConverter;

public class NamedArgumentConverter implements ArgumentConverter {
    @Override
    public Object convert(Object o, ParameterContext parameterContext) throws ArgumentConversionException {
        if (! (o instanceof Named)) {
            throw new ArgumentConversionException("Can only convert arguments of type " + Named.class.getTypeName());
        }
        return ((Named<?>) o).getPayload();
    }
}
