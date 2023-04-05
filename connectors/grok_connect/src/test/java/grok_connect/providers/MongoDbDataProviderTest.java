package grok_connect.providers;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
import com.mongodb.MongoClient;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoDatabase;
import com.mongodb.client.model.InsertOneModel;
import grok_connect.connectors_info.Credentials;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DataProvider;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameComparator;
import grok_connect.providers.utils.NamedArgumentConverter;
import grok_connect.providers.utils.Provider;
import grok_connect.utils.ProviderManager;
import grok_connect.utils.QueryMonitor;
import grok_connect.utils.SettingsManager;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.junit.jupiter.api.*;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.converter.ConvertWith;
import org.junit.jupiter.params.provider.MethodSource;
import org.mockito.Mockito;
import org.testcontainers.containers.MongoDBContainer;
import org.testcontainers.containers.wait.strategy.Wait;
import org.testcontainers.junit.jupiter.Container;
import org.testcontainers.junit.jupiter.Testcontainers;
import org.testcontainers.utility.DockerImageName;
import org.testcontainers.utility.MountableFile;
import serialization.DataFrame;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Type;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.Duration;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

@Disabled("In progress")
@TestInstance(TestInstance.Lifecycle.PER_CLASS)
@Testcontainers
class MongoDbDataProviderTest {
    private static final Provider type = Provider.MONGO_DB;
    @Container
    private static final MongoDBContainer container = new MongoDBContainer(DockerImageName.parse("mongo:5"))
            .withCopyToContainer(MountableFile.forClasspathResource("scripts/mongodb/one_line.json"),
                    "/etc/one_line.json")
            .withCopyToContainer(MountableFile.forClasspathResource("scripts/mongodb/mocks.json"),
                    "/etc/mocks.json");
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
        provider = spy.getByName(type.getProperties().get("providerName").toString());
//        MongoClient client = new MongoClient(container.getHost(), container.getFirstMappedPort());
//        MongoDatabase test = client.getDatabase("test");
//        MongoCollection<Document> mocks = test.getCollection("mocks");
//        List<InsertOneModel<Document>> docs = new ArrayList<>();
//        try {
//            List<String> list = Files.readAllLines(Paths.get("src/test/resources/scripts/mongodb/mocks.json"));
//            String join = String.join("", list);
//            String substring = join.substring(1, join.length() - 1);
//        } catch (IOException e) {
//            throw new RuntimeException(e);
//        }
    }

    @BeforeEach
    public void beforeEach() {
        Credentials credentials = new Credentials();
        credentials.parameters.put(DbCredentials.LOGIN, type.getProperties().get("user"));
        credentials.parameters.put(DbCredentials.PASSWORD, type.getProperties().get("password"));
        connection = new DataConnection();
        connection.credentials = credentials;
        connection.dataSource = provider.descriptor.type;
        connection.parameters.put(DbCredentials.SERVER, type.getProperties().get("server"));
        connection.parameters.put(DbCredentials.PORT, Double.valueOf(type.getProperties().get("port").toString()));
        connection.parameters.put(DbCredentials.DB, type.getProperties().get("db"));
    }

    @DisplayName("Tests of testConnection(DataConnection conn)")
    @Test
    public void testConnection() {
        String expected = DataProvider.CONN_AVAILABLE;
        String actual = Assertions.assertDoesNotThrow(() -> provider.testConnection(connection));
        Assertions.assertEquals(expected, actual);
    }

    @DisplayName("Mongo all types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MongoDbObjectsMother#checkOutputAllTypes_ok")
    public void checkOutputAllTypes_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Mongo string return")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MongoDbObjectsMother#checkStringReturn_ok")
    public void checkStringReturn_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Mongo one document return")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MongoDbObjectsMother#checkOneDocument_ok")
    public void checkOneDocument_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Mongo null return")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MongoDbObjectsMother#checkNullResult_ok")
    public void checkNullResult_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }
}
