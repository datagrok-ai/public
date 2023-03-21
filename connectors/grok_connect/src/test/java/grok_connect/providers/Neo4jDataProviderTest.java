package grok_connect.providers;

import grok_connect.connectors_info.Credentials;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.providers.utils.DataFrameComparator;
import grok_connect.providers.utils.Provider;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryMonitor;
import grok_connect.utils.SettingsManager;
import org.apache.log4j.Logger;
import org.junit.jupiter.api.*;
import org.mockito.Mockito;
import org.testcontainers.containers.Neo4jContainer;
import org.testcontainers.junit.jupiter.Container;
import org.testcontainers.junit.jupiter.Testcontainers;
import org.testcontainers.utility.MountableFile;
import java.io.IOException;

@Testcontainers
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class Neo4jDataProviderTest {
    private static final Provider type = Provider.NEO4J;
    private static final String HOST_SCRIPT_PATH = "scripts/neo4j/neo4j_all.cypher";
    private static final String CONTAINER_SCRIPT_PATH = "/var/lib/neo4j/neo4j_all.cypher";
    @Container
    private static final Neo4jContainer<?> container = new Neo4jContainer<>(type
            .getProperties().get("image").toString())
            .withoutAuthentication()
            .withCopyFileToContainer(MountableFile
                    .forClasspathResource(HOST_SCRIPT_PATH),
                    CONTAINER_SCRIPT_PATH);
    private JdbcDataProvider provider;
    private DataConnection connection;
    private DataFrameComparator dataFrameComparator;

    @BeforeAll
    public void init() {
        try {
            container.execInContainer("cypher-shell", "-f", CONTAINER_SCRIPT_PATH);
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException("Something went wrong when executing init script", e);
        }
        dataFrameComparator = new DataFrameComparator();
        SettingsManager settingsManager = SettingsManager.getInstance();
        settingsManager.initSettingsWithDefaults();
        Logger mockLogger = Mockito.mock(Logger.class);
        QueryMonitor mockMonitor = Mockito.mock(QueryMonitor.class);
        ProviderManager providerManager = new ProviderManager(mockLogger);
        ProviderManager spy = Mockito.spy(providerManager);
        Mockito.when(spy.getQueryMonitor()).thenReturn(mockMonitor);
        provider = spy.getByName(type.getProperties().get("providerName").toString());
    }

    @BeforeEach
    public void beforeEach() {
        Credentials credentials = new Credentials();
        connection = new DataConnection();
        connection.credentials = credentials;
        connection.dataSource = provider.descriptor.type;
        connection.parameters.put(DbCredentials.SERVER, container.getHost());
        connection.parameters.put(DbCredentials.PORT, (double) container.getFirstMappedPort());
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
}
