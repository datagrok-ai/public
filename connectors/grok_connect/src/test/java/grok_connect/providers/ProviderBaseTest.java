package grok_connect.providers;

import grok_connect.connectors_info.Credentials;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataProvider;
import grok_connect.providers.utils.Providers;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryMonitor;
import grok_connect.utils.SettingsManager;
import org.apache.log4j.Logger;
import org.junit.jupiter.api.*;
import org.mockito.Mockito;
import org.testcontainers.containers.JdbcDatabaseContainer;

/**
 * Base test class for providers, has common tests, fields and lifecycle
  */
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
public abstract class ProviderBaseTest {
    private final Providers type;
    protected JdbcDatabaseContainer<?> container;
    protected JdbcDataProvider provider;
    protected Credentials credentials;
    protected DataConnection connection;

    protected ProviderBaseTest(Providers type) {
        this.container = type.getContainer();
        this.type = type;
    }

    @BeforeAll
    public void init() {
        credentials = new Credentials();
        credentials.parameters.put("login", container.getUsername());
        credentials.parameters.put("password", container.getPassword());
        SettingsManager settingsManager = SettingsManager.getInstance();
        settingsManager.initSettingsWithDefaults();
        Logger mockLogger = Mockito.mock(Logger.class);
        QueryMonitor mockMonitor = Mockito.mock(QueryMonitor.class);
        ProviderManager providerManager = new ProviderManager(mockLogger);
        ProviderManager spy = Mockito.spy(providerManager);
        Mockito.when(spy.getQueryMonitor()).thenReturn(mockMonitor);
        provider = (JdbcDataProvider) spy.getByName(type.getProviderName());
    }

    @BeforeEach
    public void beforeEach() {
        connection = new DataConnection();
        credentials.parameters.put("login", container.getUsername());
        credentials.parameters.put("password", container.getPassword());
        connection.credentials = credentials;
        connection.dataSource = provider.descriptor.type;
        connection.parameters.put("server", container.getHost());
        connection.parameters.put("port", (double) container.getFirstMappedPort());
        connection.parameters.put("db", container.getDatabaseName());
    }

    @DisplayName("Test whether container with db is running")
    @Order(1)
    @Test
    public void database_isRunning_ok() {
        Assertions.assertTrue(container.isRunning());
    }

    @DisplayName("Tests of testConnection(DataConnection conn)")
    @Test
    public void testConnection() {
        String expected = DataProvider.CONN_AVAILABLE;
        String actual = Assertions.assertDoesNotThrow(() -> provider.testConnection(connection));
        Assertions.assertEquals(expected, actual);
    }

    // I think, it would be better to implement some Global exception handler
    // and not return String, when something went wrong in this method
    @DisplayName("Test of testConnection(DataConnection conn) when no credentials provided")
    @Test
    public void testConnection_notOk() {
        connection.credentials = null;
        String result = Assertions.assertDoesNotThrow(() -> provider.testConnection(connection));
        Assertions.assertTrue(result.startsWith("ERROR"));
    }

    @AfterAll
    public void destroy() {
        container.stop();
    }
}
