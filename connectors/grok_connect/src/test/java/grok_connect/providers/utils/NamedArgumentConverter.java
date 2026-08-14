package grok_connect.providers.utils;

import org.junit.jupiter.api.Named;
import org.junit.jupiter.api.extension.ParameterContext;
import org.junit.jupiter.params.converter.ArgumentConversionException;
import org.junit.jupiter.params.converter.ArgumentConverter;

public class NamedArgumentConverter implements ArgumentConverter {
    @Override
    public Object convert(Object o, ParameterContext parameterContext) throws ArgumentConversionException {
        // junit >= 5.11 unwraps Named before invoking converters; older versions pass it wrapped
        return o instanceof Named ? ((Named<?>) o).getPayload() : o;
    }
}
