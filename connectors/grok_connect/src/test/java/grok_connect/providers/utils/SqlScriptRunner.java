package grok_connect.providers.utils;

import org.junit.jupiter.api.extension.ExtensionContext;
import org.junit.jupiter.api.extension.InvocationInterceptor;
import org.junit.jupiter.api.extension.ReflectiveInvocationContext;
import org.testcontainers.containers.JdbcDatabaseContainer;
import org.testcontainers.ext.ScriptUtils;
import org.testcontainers.jdbc.JdbcDatabaseDelegate;
import java.lang.reflect.Field;
import java.lang.reflect.Method;

/**
 * Interceptor of test methods for running sql scripts
 */
public class SqlScriptRunner implements InvocationInterceptor {
    @Override
    public void interceptTestMethod(Invocation<Void> invocation,
                                    ReflectiveInvocationContext<Method> invocationContext,
                                    ExtensionContext extensionContext) throws Throwable {
        proceed(invocation, invocationContext,extensionContext);
    }

    @Override
    public void interceptTestTemplateMethod(Invocation<Void> invocation,
                                            ReflectiveInvocationContext<Method> invocationContext,
                                            ExtensionContext extensionContext) throws Throwable {
        proceed(invocation, invocationContext,extensionContext);
    }

    private void proceed(Invocation<Void> invocation,
                         ReflectiveInvocationContext<Method> invocationContext,
                         ExtensionContext extensionContext) throws Throwable {
        Method executable = invocationContext.getExecutable();
        if (executable.isAnnotationPresent(Sql.class)) {
            Sql annotation = executable.getAnnotation(Sql.class);
            JdbcDatabaseContainer<?> container = getContainer(extensionContext);
            JdbcDatabaseDelegate delegate = new JdbcDatabaseDelegate(container, "");
            ScriptUtils.runInitScript(delegate, annotation.restorePath());
            ScriptUtils.runInitScript(delegate, annotation.path());
            invocation.proceed();
        } else {
            invocation.proceed();
        }
    }

    private JdbcDatabaseContainer<?> getContainer(ExtensionContext extensionContext) {
        Class<?> clazz = extensionContext.getRequiredTestClass().getSuperclass();
        Object instance = extensionContext.getRequiredTestInstance();
        return findField(clazz, instance);
    }

    private JdbcDatabaseContainer<?> findField(Class<?> clazz, Object instance) {
        while(clazz != null) {
            Field[] declaredFields = clazz.getDeclaredFields();
            for (Field field : declaredFields) {
                if (field.getType().isAssignableFrom(JdbcDatabaseContainer.class)) {
                    try {
                        field.setAccessible(true);
                        return (JdbcDatabaseContainer<?>) field.get(instance);
                    } catch (IllegalAccessException e) {
                        throw new RuntimeException("Something went wrong "
                                + "when getting field of type <JdbcDatabaseContainer.class>", e);
                    }
                }
            }
            clazz = clazz.getSuperclass();
        }
        throw new UnsupportedOperationException("Class declaring <Sql> annotation over method "
                + "should have field of type <JdbcDatabaseContainer.class>");
    }
}
