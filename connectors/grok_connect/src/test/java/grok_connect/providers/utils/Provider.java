package grok_connect.providers.utils;

import org.testcontainers.containers.*;
import org.testcontainers.containers.wait.strategy.DockerHealthcheckWaitStrategy;
import org.testcontainers.utility.DockerImageName;
import java.io.FileReader;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Properties;

/**
 * Enum that contains all necessary data related to specific provider and it's container
 */
public enum Provider {
    PI("src/test/resources/properties/pi.properties"),

    MONGO_DB("src/test/resources/properties/mongodb.properties"),

    CLICKHOUSE("src/test/resources/properties/clickhouse.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            container = new ClickHouseContainer(DockerImageName.parse(properties.get("image").toString()))
                    .waitingFor(new DockerHealthcheckWaitStrategy())
                    .withInitScript(properties.get("initScript").toString());
            container.start();
            return container;
        }
    },

    MARIADB("src/test/resources/properties/mariadb.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            container = new MariaDBContainer<>()
                    .withDatabaseName(properties.get("database").toString())
                    .withUsername(properties.get("user").toString())
                    .withPassword(properties.get("password").toString())
                    .withInitScript(properties.get("initScript").toString());
            container.start();
            return container;
        }
    },

    DB2("src/test/resources/properties/db2.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {

            container = new Db2Container(properties.get("image").toString())
                    .withDatabaseName(properties.get("database").toString())
                    .withUsername(properties.get("user").toString())
                    .withPassword(properties.get("password").toString())
                    .withInitScript(properties.get("initScript").toString())
                    .acceptLicense();
            container.start();
            return container;
        }
    },

    NEO4J("src/test/resources/properties/neo4j.properties"),

    HIVE2("src/test/resources/properties/hive2.properties"),

    REDSHIFT("src/test/resources/properties/redshift.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {

            container = new MySQLContainer<>(properties.get("image").toString())
                    .withDatabaseName(properties.get("database").toString())
                    .withUsername(properties.get("user").toString())
                    .withPassword(properties.get("password").toString())
                    .withEnv("MYSQL_ROOT_PASSWORD", "datagrok")
                    .withInitScript(properties.get("initScript").toString());
            container.start();
            return container;
        }
    },

    MYSQL("src/test/resources/properties/mysql.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            container = new MySQLContainer<>(properties.get("image").toString())
                    .withDatabaseName(properties.get("database").toString())
                    .withUsername(properties.get("user").toString())
                    .withPassword(properties.get("password").toString())
                    .withEnv("MYSQL_ROOT_PASSWORD", "datagrok")
                    .withInitScript(properties.get("initScript").toString());
            container.start();
            return container;
        }
    },

    ATHENA("src/test/resources/properties/athena.properties"),

    MSSQL("src/test/resources/properties/mssql.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            container = new MSSQLServerContainer(properties.get("image").toString())
                    .withInitScript(properties.get("initScript").toString());
            container.start();
            return container;
        }
    },

    SNOWFLAKE("src/test/resources/properties/snowflake.properties"),

    ORACLE("src/test/resources/properties/oracle.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            container = new OracleContainer(properties.get("image").toString())
                    .withDatabaseName(properties.get("database").toString())
                    .withUsername(properties.get("user").toString())
                    .withPassword(properties.get("password").toString())
                    .withInitScript(properties.get("initScript").toString());
            container.start();
            return container;
        }
    },

    POSTGRESQL("src/test/resources/properties/postgresql.properties") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            container = new PostgreSQLContainer<>(properties.get("image").toString())
                    .withDatabaseName(properties.get("database").toString())
                    .withUsername(properties.get("user").toString())
                    .withPassword(properties.get("password").toString())
                    .withInitScript(properties.get("initScript").toString())
                    .withClasspathResourceMapping(properties.get("volumePath").toString(),
                                "/etc/", BindMode.READ_ONLY);
            container.start();
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

    public static Provider fromName(String providerName) {
        for (Provider provider: Provider.values()) {
            String currentProviderName = provider.getProperties()
                    .get("providerName")
                    .toString().replaceAll(" ", "");
            if (currentProviderName.equalsIgnoreCase(providerName)) {
                return provider;
            }
        }
        throw new NoSuchElementException("No such provider is registered");
    }

    public JdbcDatabaseContainer<?> getContainer() {
        if (this.container == null) {
            return newJdbcContainer();
        }
        return container;
    }

    public Properties getProperties() {
        return properties;
    }

    protected JdbcDatabaseContainer<?> newJdbcContainer() {
        throw new UnsupportedOperationException("Override method newJdbcContainer()");
    }
}
