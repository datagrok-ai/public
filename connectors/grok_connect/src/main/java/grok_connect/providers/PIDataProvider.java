package grok_connect.providers;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncParam;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import grok_connect.utils.ProviderManager;

public class PIDataProvider extends JdbcDataProvider {
    public PIDataProvider(ProviderManager providerManager) {
        super(providerManager);
        driverClassName = "com.osisoft.jdbc.Driver";
        descriptor = new DataSource();
        descriptor.type = "PI";
        descriptor.description = "Query PI database";
        descriptor.connectionTemplate = new ArrayList<Property>() {{
            add(new Property(Property.STRING_TYPE, DbCredentials.ACCESS_SERVER));
            add(new Property(Property.STRING_TYPE, DbCredentials.SERVER));
            add(new Property(Property.STRING_TYPE, DbCredentials.INITIAL_CATALOG, "Name of the catalog,"
                    + " e.g. piarchive, pibatch, pids, pifunction, piheading, "
                    + "pilog, pimodule, pipoint, pisystem or piuser"));
            add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                    DbCredentials.CONNECTION_STRING_DESCRIPTION, new Prop("textarea")));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_SCHEMA));
            add(new Property(Property.BOOL_TYPE, DbCredentials.CACHE_RESULTS));
            add(new Property(Property.STRING_TYPE, DbCredentials.CACHE_INVALIDATE_SCHEDULE));
        }};
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
        descriptor.nameBrackets = "\"";
        descriptor.defaultSchema = "public";
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        return "jdbc:pioledb://" + conn.parameters.get(DbCredentials.ACCESS_SERVER) + "/Data Source="
                + conn.getServer() + ";Initial Catalog=" + conn.get(DbCredentials.INITIAL_CATALOG)
                + ";Integrated Security=SSPI";
    }

    @Override
    protected void appendQueryParam(DataQuery dataQuery, String paramName, StringBuilder queryBuffer) {
        FuncParam param = dataQuery.getParam(paramName);
        if (param.propertyType.equals("list")) {
            @SuppressWarnings("unchecked")
            List<String> values = ((ArrayList<String>) param.value);
            queryBuffer.append(values.stream().map(value -> "?").collect(Collectors.joining(", ")));
        } else {
            queryBuffer.append("?");
        }
    }

    @Override
    protected int setArrayParamValue(PreparedStatement statement, int n, FuncParam param) throws SQLException {
        @SuppressWarnings (value="unchecked")
        ArrayList<Object> lst = (ArrayList<Object>)param.value;
        if (lst == null || lst.size() == 0) {
            statement.setObject(n, null);
            return 0;
        }
        for (int i = 0; i < lst.size(); i++) {
            statement.setObject(n + i, lst.get(i));
        }
        return lst.size() - 1;
    }
}
