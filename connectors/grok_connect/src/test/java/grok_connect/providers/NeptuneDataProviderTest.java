package grok_connect.providers;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.*;
import grok_connect.providers.utils.DataFrameComparator;
import grok_connect.providers.utils.NamedArgumentConverter;
import grok_connect.providers.utils.Provider;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryMonitor;
import grok_connect.utils.SettingsManager;
import org.junit.jupiter.api.*;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.converter.ConvertWith;
import org.junit.jupiter.params.provider.MethodSource;
import org.mockito.Mockito;
import serialization.DataFrame;

@Disabled("Until test instance of Neptune will be available")
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
class NeptuneDataProviderTest {
    private static final Provider type = Provider.NEPTUNE;
    private JdbcDataProvider provider;
    private DataConnection connection;
    private DataFrameComparator dataFrameComparator;

    @BeforeAll
    public void init() {
        dataFrameComparator = new DataFrameComparator();
        SettingsManager settingsManager = SettingsManager.getInstance();
        settingsManager.initSettingsWithDefaults();
        QueryMonitor mockMonitor = Mockito.mock(QueryMonitor.class);
        ProviderManager providerManager = new ProviderManager();
        ProviderManager spy = Mockito.spy(providerManager);
        Mockito.when(spy.getQueryMonitor()).thenReturn(mockMonitor);
        GrokConnect.providerManager = providerManager;
        provider = spy.getByName(type.getProperties().get("providerName").toString());
    }

    @BeforeEach
    public void beforeEach() {
        Credentials credentials = new Credentials();
        credentials.parameters.put(DbCredentials.ACCESS_KEY, type.getProperties().get("accessKey"));
        credentials.parameters.put("secretAccessKey", type.getProperties().get("secretAccessKey"));
        connection = new DataConnection();
        connection.credentials = credentials;
        connection.dataSource = provider.descriptor.type;
        connection.parameters.put(DbCredentials.CONNECTION_STRING, type.getProperties().get("connectionString"));
        connection.parameters.put("serviceRegion", type.getProperties().get("serviceRegion"));
    }

    @DisplayName("Tests of testConnection(DataConnection conn)")
    @Test
    public void testConnection() {
        String expected = DataProvider.CONN_AVAILABLE;
        String actual = Assertions.assertDoesNotThrow(() -> provider.testConnection(connection));
        Assertions.assertEquals(expected, actual);
    }

    @DisplayName("Some test query")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.NeptuneObjectsMother#someTest")
    public void someTest(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }
}
