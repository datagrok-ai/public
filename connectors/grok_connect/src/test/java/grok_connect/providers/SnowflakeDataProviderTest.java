package grok_connect.providers;

import grok_connect.connectors_info.Credentials;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameComparator;
import grok_connect.providers.utils.Provider;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryMonitor;
import grok_connect.utils.SettingsManager;
import org.apache.log4j.Logger;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInstance;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.mockito.Mockito;
import serialization.DataFrame;

@TestInstance(TestInstance.Lifecycle.PER_CLASS)
class SnowflakeDataProviderTest {
    private static final Provider type = Provider.SNOWFLAKE;
    private JdbcDataProvider provider;
    private DataConnection connection;
    private DataFrameComparator dataFrameComparator;

    @BeforeAll
    public void init() {
        dataFrameComparator = new DataFrameComparator();
        SettingsManager settingsManager = SettingsManager.getInstance();
        settingsManager.initSettingsWithDefaults();
        Logger mockLogger = Mockito.mock(Logger.class);
        QueryMonitor mockMonitor = Mockito.mock(QueryMonitor.class);
        ProviderManager providerManager = new ProviderManager(mockLogger);
        ProviderManager spy = Mockito.spy(providerManager);
        Mockito.when(spy.getQueryMonitor()).thenReturn(mockMonitor);
        provider = (JdbcDataProvider) spy.getByName(type.getProperties().get("providerName").toString());
    }

    @BeforeEach
    public void beforeEach() {
        Credentials credentials = new Credentials();
        credentials.parameters.put(DbCredentials.LOGIN, type.getProperties().get("user"));
        credentials.parameters.put(DbCredentials.PASSWORD, type.getProperties().get("password"));
        connection = new DataConnection();
        connection.credentials = credentials;
        connection.dataSource = provider.descriptor.type;
        connection.parameters.put(DbCredentials.ACCOUNT_LOCATOR, type.getProperties().get("accountLocator"));
        connection.parameters.put(DbCredentials.REGION_ID, type.getProperties().get("region"));
        connection.parameters.put(DbCredentials.DB, type.getProperties().get("database"));
        connection.parameters.put(DbCredentials.CLOUD, type.getProperties().get("cloud"));
        connection.parameters.put(DbCredentials.WAREHOUSE, type.getProperties().get("warehouse"));
    }

    @DisplayName("Tests of testConnection(DataConnection conn)")
    @Test
    public void testConnection() {
        String expected = DataProvider.CONN_AVAILABLE;
        String actual = Assertions.assertDoesNotThrow(() -> provider.testConnection(connection));
        Assertions.assertEquals(expected, actual);
    }

    @DisplayName("Test of testConnection(DataConnection conn) when no credentials provided")
    @Test
    public void testConnection_notOk() {
        connection.credentials = null;
        String result = Assertions.assertDoesNotThrow(() -> provider.testConnection(connection));
        Assertions.assertTrue(result.startsWith("ERROR"));
    }

    @DisplayName("Test of getSchemas() method with correct DataConnection")
    @ParameterizedTest(name = "CORRECT ARGUMENTS")
    @MethodSource("grok_connect.providers.arguments_provider.SnowflakeObjectsMother#getSchemas_ok")
    public void getSchemas_ok(DataFrame expected) {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getSchemas(connection));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @Disabled
    @Test
    public void getSchemas_notOk() {
        // method probably should throw something when bad input
    }

    @DisplayName("Test of getSchema() method with correct DataConnection")
    @ParameterizedTest(name = "CORRECT ARGUMENTS")
    @MethodSource("grok_connect.providers.arguments_provider.SnowflakeObjectsMother#getSchema_ok")
    public void getSchema_ok(DataFrame expected) {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getSchema(connection,
                "PUBLIC", "MOCK_DATA"));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @Disabled
    @Test
    public void getSchema_notOk() {
        // method probably should throw something when bad input
    }

    @DisplayName("Parameters support")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource({"grok_connect.providers.arguments_provider.CommonObjectsMother#checkParameterSupport_ok",
            "grok_connect.providers.arguments_provider.CommonObjectsMother#checkMultipleParametersSupport_ok",
            "grok_connect.providers.arguments_provider.CommonObjectsMother#checkListParameterSupport_ok",
            "grok_connect.providers.arguments_provider.CommonObjectsMother#checkRegexSupport_ok"})
    public void checkParameterSupport_ok(FuncCall funcCall, DataFrame expected) {
        prepareDataFrame(expected);
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support for datetime")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.CommonObjectsMother#checkDatesParameterSupport_ok")
    public void checkDatesParameterSupport_ok(FuncCall funcCall, DataFrame expected) {
        prepareDataFrame(expected);
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for snowflake date, time, timestamp types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.SnowflakeObjectsMother#checkOutputDataFrame_dateTypes_ok")
    public void checkOutputDataFrame_dateTypes_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for snowflake numeric types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.SnowflakeObjectsMother#checkOutputDataFrame_numericTypes_ok")
    public void checkOutputDataFrame_numericTypes_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for snowflake binary type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.SnowflakeObjectsMother#checkOutputDataFrame_binaryType_ok")
    public void checkOutputDataFrame_binaryType_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for snowflake geo type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.SnowflakeObjectsMother#checkOutputDataFrame_geoType_ok")
    public void checkOutputDataFrame_geoType_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for snowflake semi-structured types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.SnowflakeObjectsMother#checkOutputDataFrame_semiStructuredTypes_ok")
    public void checkOutputDataFrame_semiStructuredTypes_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    private void prepareDataFrame(DataFrame dataFrame) {
        dataFrame.columns.forEach(column -> column.name = column.name.toUpperCase());
    }
}