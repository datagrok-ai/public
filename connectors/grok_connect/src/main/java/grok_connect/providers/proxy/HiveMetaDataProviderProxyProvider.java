package grok_connect.providers.proxy;

import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.ProviderManager;
import javassist.util.proxy.MethodHandler;
import javassist.util.proxy.ProxyFactory;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class HiveMetaDataProviderProxyProvider {
    private static final List<String> SUPPORTED_METHODS =
            Collections.unmodifiableList(Arrays.asList("getSchemaSql", "getSchemasSql"));

    public JdbcDataProvider getProxy(ProviderManager providerManager,
                                     Class<?> clazz) {
        ProxyFactory factory = new ProxyFactory();
        factory.setSuperclass(clazz);
        factory.setFilter((method -> SUPPORTED_METHODS.contains(method.getName())));
        try {
            return (JdbcDataProvider) factory.create(new Class<?>[]{ProviderManager.class},
                    new Object[]{providerManager}, getHandler());
        } catch (NoSuchMethodException | InstantiationException | IllegalAccessException | InvocationTargetException e) {
            throw new RuntimeException("Something went wrong when creating proxy for "
                    + factory.getSuperclass().getName(), e);
        }
    }

    private MethodHandler getHandler() {
        return (self, thisMethod, proceed, args) -> {
            if (thisMethod.getName().equals("getSchemaSql")) {
                return handleGetSchemaImpl(args);
            }
            return handleGetSchemasImpl(args);
        };
    }

    private String handleGetSchemasImpl(Object[] args) {
        return "SELECT DISTINCT \"NAME\" AS table_schema  FROM \"DBS\";";
    }

    private String handleGetSchemaImpl(Object[] args) {
        Object schema = args[1];
        Object table = args[2];
        return String.format("SELECT \"DBS\".\"NAME\" as table_schema, \"TBL_NAME\" "
                + "as table_name, \"COLUMN_NAME\" as column_name, \"TYPE_NAME\" as data_type, \n"
                + "CASE WHEN \"TBL_TYPE\" = 'VIRTUAL_VIEW' THEN 1 ELSE 0 END AS is_view FROM \"TBLS\", "
                + "\"COLUMNS_V2\", \"DBS\" WHERE \"TBL_ID\"=\"CD_ID\" AND \"DBS\".\"DB_ID\" = \"TBLS\".\"DB_ID\"" +
                "%s%s;", schema != null ? String.format(" AND \"DBS\".\"NAME\" = '%s'%s",schema, table != null ? " AND " : "") : "",
                table != null ? String.format("\"TBL_NAME\" = '%s'", table) : "");
    }
}
