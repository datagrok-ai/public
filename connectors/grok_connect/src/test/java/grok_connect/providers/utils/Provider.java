package grok_connect.providers.utils;

import org.testcontainers.containers.BindMode;
import org.testcontainers.containers.JdbcDatabaseContainer;
import org.testcontainers.containers.MSSQLServerContainer;
import org.testcontainers.containers.OracleContainer;
import org.testcontainers.containers.PostgreSQLContainer;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;

/**
 * Enum that contains all necessary data related to specific provider and it's container
 */
public enum Provider {
    ATHENA("src/test/resources/properties/athena.properties"),
    MSSQL("src/test/resources/properties/mssql.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            if (this.container == null) {
                container = new MSSQLServerContainer(properties.get("image").toString())
                        .withInitScript(properties.get("initScript").toString());
                container.start();
            }
            return container;
        }
    },

    SNOWFLAKE("src/test/resources/properties/snowflake.properties"),

    ORACLE("src/test/resources/properties/oracle.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            if (this.container == null) {
                container = new OracleContainer(properties.get("image").toString())
                        .withDatabaseName(properties.get("database").toString())
                        .withUsername(properties.get("user").toString())
                        .withPassword(properties.get("password").toString())
                        .withInitScript(properties.get("initScript").toString());
                container.start();
            }
            return container;
        }
    },

    POSTGRESQL("src/test/resources/properties/postgresql.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            if (this.container == null) {
                container = new PostgreSQLContainer<>(properties.get("image").toString())
                        .withDatabaseName(properties.get("database").toString())
                        .withUsername(properties.get("user").toString())
                        .withPassword(properties.get("password").toString())
                        .withClasspathResourceMapping(properties.get("volumePath").toString(),
                                "/etc/", BindMode.READ_ONLY);
                container.start();
            }
            return container;
        }
    };

    protected final Properties properties;
    protected JdbcDatabaseContainer<?> container;

    Provider(String propertyPath) {
        Properties properties = new Properties();
        try {
            properties.load(new FileReader(propertyPath));
        } catch (IOException e) {
            throw new RuntimeException("Something went wrong when reading from " + propertyPath, e);
        }
        this.properties = properties;
    }

    protected JdbcDatabaseContainer<?> newJdbcContainer() {
        throw new UnsupportedOperationException("Override method newJdbcContainer()");
    }

    public JdbcDatabaseContainer<?> getContainer() {
        return newJdbcContainer();
    }

    public Properties getProperties() {
        return properties;
    }
}
