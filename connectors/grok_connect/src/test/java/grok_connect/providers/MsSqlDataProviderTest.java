package grok_connect.providers;

import grok_connect.connectors_info.Credentials;
import grok_connect.connectors_info.DataConnection;
import grok_connect.connectors_info.DbCredentials;
import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.NamedArgumentConverter;
import grok_connect.providers.utils.Provider;
import grok_connect.providers.utils.Sql;
import grok_connect.providers.utils.SqlScriptRunner;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.converter.ConvertWith;
import org.junit.jupiter.params.provider.MethodSource;
import serialization.DataFrame;

@ExtendWith(SqlScriptRunner.class)
class MsSqlDataProviderTest extends ContainerizedProviderBaseTest {
    private static final String DEFAULT_DATABASE_NAME = "master";
    private static final String DEFAULT_SCHEMA_NAME = "dbo";

    protected MsSqlDataProviderTest(Provider provider) {
        super(provider);
    }

    @Override
    @BeforeEach
    public void beforeEach() {
        credentials = new Credentials();
        credentials.parameters.put(DbCredentials.LOGIN, container.getUsername());
        credentials.parameters.put(DbCredentials.PASSWORD, container.getPassword());
        connection = new DataConnection();
        connection.credentials = credentials;
        connection.dataSource = provider.descriptor.type;
        connection.parameters.put(DbCredentials.SERVER, container.getHost());
        connection.parameters.put(DbCredentials.PORT, (double) container.getFirstMappedPort());
        connection.parameters.put(DbCredentials.DB, DEFAULT_DATABASE_NAME); // CAN'T GET DB NAME FROM MSSQL CONTAINER
        System.out.println(container.getJdbcUrl());
    }

    @DisplayName("Test of getCatalogs() method with correct DataConnection")
    @Test
    public void getCatalogs_ok() {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getCatalogs(connection));
        Assertions.assertNotNull(actual);
        Assertions.assertEquals(1, actual.columns.size());
        Assertions.assertEquals("catalog_name", actual.columns.get(0).name);
        Assertions.assertTrue(actual.rowCount > 0);
    }

    @DisplayName("Test of getSchemas() method with correct DataConnection")
    @ParameterizedTest(name = "CORRECT ARGUMENTS")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#getSchemas_ok")
    @Sql(path = "scripts/mssql/mssql_basic_types.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void getSchemas_ok(DataFrame expected) {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getSchemas(connection));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @Disabled
    @Test
    public void getSchemas_notOk() {
        // method probably should throw something when bad input
    }

    @Disabled
    @DisplayName("Test of getSchema() method with correct DataConnection")
    @ParameterizedTest(name = "CORRECT ARGUMENTS")
    @Sql(path = "scripts/mssql/mssql_basic_types.sql",
            restorePath = "scripts/mssql/drop.sql")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#getSchema_ok")
    public void getSchema_ok(DataFrame expected) {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getSchema(connection,
                DEFAULT_SCHEMA_NAME, "MOCK_DATA", false));
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
            "grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkMultipleParametersSupport_ok",
            "grok_connect.providers.arguments_provider.CommonObjectsMother#checkListParameterSupport_ok"})
    @Sql(path = "scripts/mssql/mssql_basic_types.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkParameterSupport_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        prepareDataFrame(expected);
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support for datetime")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.CommonObjectsMother#checkDatesParameterSupport_ok")
    @Sql(path = "scripts/mssql/mssql_dates_patterns.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkDatesParameterSupport_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for SQL SERVER date types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkOutputDataFrame_dateTypes_ok")
    @Sql(path = "scripts/mssql/mssql_date_types.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkOutputDataFrame_dateTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for SQL SERVER character types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkOutputDataFrame_characterTypes_ok")
    @Sql(path = "scripts/mssql/mssql_character_types.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkOutputDataFrame_characterTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for SQL SERVER xml type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkOutputDataFrame_xmlType_ok")
    @Sql(path = "scripts/mssql/mssql_xml_type.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkOutputDataFrame_xmlType_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for SQL SERVER binary types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkOutputDataFrame_binaryTypes_ok")
    @Sql(path = "scripts/mssql/mssql_binary_types.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkOutputDataFrame_binaryTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for SQL SERVER numeric types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkOutputDataFrame_numericTypes_ok")
    @Sql(path = "scripts/mssql/mssql_numeric_types.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkOutputDataFrame_numericTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for SQL SERVER money types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkOutputDataFrame_moneyTypes_ok")
    @Sql(path = "scripts/mssql/mssql_money_type.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkOutputDataFrame_moneyTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for SQL SERVER spatial types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.MsSqlObjectsMother#checkOutputDataFrame_spatialTypes_ok")
    @Sql(path = "scripts/mssql/mssql_spatial_identifier.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkOutputDataFrame_spatialTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("SQL SERVER null safety")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.CommonObjectsMother#checkNullSupport_ok")
    @Sql(path = "scripts/mssql/mssql_null.sql",
            restorePath = "scripts/mssql/drop.sql")
    public void checkNullSupport_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall) {
        funcCall.func.connection = connection;
        Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
    }

    private void prepareDataFrame(DataFrame dataFrame) {
        // in order to save time reuse some common's
        dataFrame.columns.removeIf(column -> column.name.equals("bool")); // mssql doesn't have boolean type
    }
}
