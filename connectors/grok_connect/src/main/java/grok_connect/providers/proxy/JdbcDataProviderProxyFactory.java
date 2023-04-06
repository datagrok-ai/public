package grok_connect.providers.proxy;

import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.ProviderManager;
import javassist.util.proxy.MethodHandler;
import javassist.util.proxy.ProxyFactory;
import org.slf4j.Logger;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;

public class JdbcDataProviderProxyFactory {
    private static final String LOGGER_FIELD_NAME = "logger";
    public JdbcDataProvider getProxy(ProviderManager providerManager,
                                     Class<?> clazz) {
        ProxyFactory factory = new ProxyFactory();
        factory.setSuperclass(clazz);
        factory.setFilter(method -> true); // log all methods
        try {
            return (JdbcDataProvider) factory.create(new Class<?>[]{ProviderManager.class},
                    new Object[]{providerManager}, getHandler());
        } catch (NoSuchMethodException | InstantiationException | IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException("Something went wrong when creating proxy for %s"
                    + factory.getSuperclass().getName(), e);
        }
    }

    private MethodHandler getHandler() {
        return (self, thisMethod, proceed, args) -> {
            Field loggerField = self.getClass().getSuperclass().getSuperclass().getDeclaredField("logger");
            loggerField.setAccessible(true);
            Logger logger = (Logger) loggerField.get(self);
            loggerField.setAccessible(false);
            logger.warn("{} was called with args {}", thisMethod.getName(), Arrays.toString(args));
            Object result;
            try {
                result = thisMethod.invoke(self, args);
            } catch (IllegalAccessException | InvocationTargetException e) {
                throw new RuntimeException(String.format("Something went wrong when invoking method %s", thisMethod), e);
            }
            logger.warn("{} was successfully proceeded", thisMethod.getName());
            return result;
        };
    }
}
