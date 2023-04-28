package grok_connect.providers.utils;

import org.junit.jupiter.api.extension.ExtensionContext;
import org.junit.jupiter.api.extension.ParameterContext;
import org.junit.jupiter.api.extension.ParameterResolutionException;
import org.junit.jupiter.api.extension.ParameterResolver;

public class ConstructorParameterResolver implements ParameterResolver {
    private static final String SUFFIX = "DataProviderTest";

    @Override
    public boolean supportsParameter(ParameterContext parameterContext, ExtensionContext extensionContext)
            throws ParameterResolutionException {
        return parameterContext.getParameter().getType() == Provider.class;
    }

    @Override
    public Object resolveParameter(ParameterContext parameterContext, ExtensionContext extensionContext)
            throws ParameterResolutionException {
        String declaringClassName = parameterContext.getDeclaringExecutable()
                .getDeclaringClass().getSimpleName();
        String dataProviderName = declaringClassName.replace(SUFFIX, "");
        return Provider.fromName(dataProviderName);
    }
}
