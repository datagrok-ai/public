package grok_connect.providers;

import grok_connect.connectors_info.Credentials;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DbCredentials;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class DatabricksProviderTest {
    @Test
    void schemaGoesIntoUrlPath() {
        DatabricksProvider provider = new DatabricksProvider();
        DataConnection conn = new DataConnection();
        conn.credentials = new Credentials();
        conn.parameters.put("workspaceURL", "dbc-1.cloud.databricks.com");
        Assertions.assertTrue(provider.getConnectionStringImpl(conn).endsWith(":443/default"));
        conn.parameters.put(DbCredentials.SCHEMA, "sales");
        Assertions.assertTrue(provider.getConnectionStringImpl(conn).endsWith(":443/sales"));
        Assertions.assertNull(provider.getProperties(conn).getProperty("ConnSchema"));
    }
}
