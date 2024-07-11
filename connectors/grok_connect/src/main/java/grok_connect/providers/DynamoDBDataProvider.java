package grok_connect.providers;

import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataSource;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.utils.Prop;
import grok_connect.utils.Property;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Properties;

public class DynamoDBDataProvider extends JdbcDataProvider {
    private static final String DEFAULT_DOMAIN = "amazonaws.com";
    private static final List<String> AVAILABLE_REGIONS = Collections.unmodifiableList(Arrays.asList(
            "OHIO", "NORTHERNVIRGINIA", "NORTHERNCALIFORNIA", "OREGON", "CAPETOWN", "HONGKONG", "JAKARTA", "MUMBAI",
            "OSAKA", "SEOUL", "SINGAPORE", "SYDNEY", "TOKYO", "CENTRAL", "BEIJING", "NINGXIA", "FRANKFURT", "IRELAND",
            "LONDON", "MILAN", "PARIS", "STOCKHOLM", "ZURICH", "BAHRAIN", "UAE", "SAOPAULO", "GOVCLOUDEAST",
            "GOVCLOUDWEST"
    ));

    public DynamoDBDataProvider() {
        driverClassName = "cdata.jdbc.amazondynamodb.AmazonDynamoDBDriver";

        descriptor = new DataSource();
        descriptor.type = "DynamoDB";
        descriptor.description = "Query DynamoDB database";
        descriptor.connectionTemplate = new ArrayList<>();
        descriptor.connectionTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.DOMAIN,
                "If not specified default domain will be used - 'amazonaws.com'"));
        Property regions = new Property(Property.STRING_TYPE, DbCredentials.REGION_ID, "AWS Region where your "
                + "DynamoDB instance is located.");
        regions.choices = new ArrayList<>();
        regions.choices.addAll(AVAILABLE_REGIONS);
        descriptor.connectionTemplate.add(regions);
        descriptor.connectionTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.CONNECTION_STRING,
                "When specified, this connection string overrides all other parameters. Example: "
                        + "'jdbc:amazondynamodb:AWSAccessKey=[Your access key];AWSSecretKey=[Your secret key];"
                        + "Domain=[Your domain];AWSRegion=[Your region];'" , new Prop("textarea")));
        descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
        descriptor.credentialsTemplate = new ArrayList<>();
        descriptor.credentialsTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.ACCESS_KEY,
                    "Set this property to the access\n"
                            + "key provided by your AWS account.\n"));
        descriptor.credentialsTemplate.add(new Property(Property.STRING_TYPE, DbCredentials.SECRET_KEY,
                    "Set this property to the secret\n"
                            + "key provided by your AWS account.", new Prop("password")));
            descriptor.nameBrackets = "\"";
    }

    @Override
    public Properties getProperties(DataConnection conn) {
        java.util.Properties properties = defaultConnectionProperties(conn);
        if (!conn.hasCustomConnectionString() && conn.ssl()) {
            properties.setProperty("ssl", "true");
        }
        return properties;
    }

    @Override
    public String getConnectionStringImpl(DataConnection conn) {
        String domain = conn.get(DbCredentials.DOMAIN);
        String region = conn.get(DbCredentials.REGION_ID);
        String accessKey = (String) conn.credentials.parameters.get(DbCredentials.ACCESS_KEY);
        String secretKey = (String) conn.credentials.parameters.get(DbCredentials.SECRET_KEY);
        return String.format("jdbc:amazondynamodb:AWSAccessKey=%s;"
                + "AWSSecretKey=%s;Domain=%s;AWSRegion=%s;", accessKey, secretKey,
                domain == null || domain.isEmpty() ? DEFAULT_DOMAIN : domain,
                region);
    }
}
