package grok_connect.table_query;

import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.GrokConnectUtil;

/** Shared identifier-quoting helpers used by both TableQuery and table_mutation SQL emission. */
public class SqlNames {
    public static String fullTableName(String tableName, String schema, String catalog, JdbcDataProvider provider) {
        if (tableName == null || provider == null)
            throw new IllegalArgumentException("tableName and provider must not be null");

        int dotCount = tableName.length() - tableName.replace(".", "").length();

        // Already fully qualified (catalog.schema.table)
        if (dotCount >= 2)
            return provider.addBrackets(tableName);

        boolean hasSchema = GrokConnectUtil.isNotEmpty(schema);
        boolean hasCatalog = GrokConnectUtil.isNotEmpty(catalog);

        // Already has schema.table
        if (dotCount == 1) {
            if (provider.descriptor.supportCatalogs && hasCatalog)
                return provider.addBrackets(catalog + "." + tableName);
            return provider.addBrackets(tableName);
        }

        // Simple table name (no dots)
        if (provider.descriptor.supportCatalogs && hasCatalog && hasSchema)
            return provider.addBrackets(catalog + "." + schema + "." + tableName);
        if (hasSchema)
            return provider.addBrackets(schema + "." + tableName);
        return provider.addBrackets(tableName);
    }
}
