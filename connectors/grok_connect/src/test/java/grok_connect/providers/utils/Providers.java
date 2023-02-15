package grok_connect.providers.utils;

import org.testcontainers.containers.BindMode;
import org.testcontainers.containers.JdbcDatabaseContainer;
import org.testcontainers.containers.PostgreSQLContainer;

/**
 * Enum that contains all necessary data related to specific provider and it's container
 */
public enum Providers {
    POSTGRESQL("Postgres", "postgres:12-alpine",
            "datagrok", "postgres", "postgres",
            "scripts/postgres/file.txt") {
        @Override
        protected JdbcDatabaseContainer<?> newJdbcContainer() {
            if (this.container == null) {
                container = new PostgreSQLContainer<>(getDefaultImage())
                        .withDatabaseName(getDatabase())
                        .withUsername(getDefaultUser())
                        .withPassword(getDefaultPassword())
                        .withClasspathResourceMapping(getVolumePath(),
                                "/etc/", BindMode.READ_ONLY);
                container.start();
            }
            return container;
        }
    };

    private final String providerName;
    private final String defaultImage;
    private final String database;
    private final String defaultUser;
    private final String defaultPassword;
    private final String volumePath;
    protected JdbcDatabaseContainer<?> container;

    Providers(String providerName, String image, String database, String defaultUser,
              String defaultPassword, String volumePath) {
        this.providerName = providerName;
        this.defaultImage = image;
        this.database = database;
        this.defaultUser = defaultUser;
        this.defaultPassword = defaultPassword;
        this.volumePath = volumePath;
    }

    protected JdbcDatabaseContainer<?> newJdbcContainer() {
        throw new UnsupportedOperationException("Override method newJdbcContainer()");
    }

    public JdbcDatabaseContainer<?> getContainer() {
        return newJdbcContainer();
    }

    public String getProviderName() {
        return providerName;
    }

    public String getDefaultImage() {
        return defaultImage;
    }

    public String getDatabase() {
        return database;
    }

    public String getDefaultUser() {
        return defaultUser;
    }

    public String getDefaultPassword() {
        return defaultPassword;
    }

    public String getVolumePath() {
        return volumePath;
    }
}
