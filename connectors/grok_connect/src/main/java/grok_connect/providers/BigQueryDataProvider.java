package grok_connect.providers;

import java.util.*;

import grok_connect.resultset.ResultSetManager;
import grok_connect.utils.*;
import grok_connect.connectors_info.*;


public class BigQueryDataProvider extends JdbcDataProvider {
    public BigQueryDataProvider(ResultSetManager resultSetManager, ProviderManager providerManager) {
        super(resultSetManager, providerManager);
        driverClassName = "net.starschema.clouddb.jdbc.BQDriver";

        descriptor = new DataSource();
        descriptor.type = "BigQuery";
        descriptor.description = "Query BigQuery database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.STRING_TYPE, "project_id", "ID of project"));
        }};
        descriptor.credentialsTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.LOGIN));
            add(new Property(Property.STRING_TYPE, DbCredentials.PASSWORD, new Prop("password")));
        }};
    }

    public String getConnectionStringImpl(DataConnection conn) {
        return "jdbc:BQDriver:" + conn.parameters.get("project_id") + "?withServiceAccount=true";
    }
}
