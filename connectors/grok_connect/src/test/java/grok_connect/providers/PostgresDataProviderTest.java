package grok_connect.providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.NamedArgumentConverter;
import grok_connect.providers.utils.Provider;
import grok_connect.providers.utils.Sql;
import grok_connect.providers.utils.SqlScriptRunner;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.converter.ConvertWith;
import org.junit.jupiter.params.provider.MethodSource;
import serialization.DataFrame;

/**
 * Test class for PostgresDataProvider
 */
@ExtendWith(SqlScriptRunner.class)
class PostgresDataProviderTest extends ContainerizedProviderBaseTest {
    private static final String INIT_SCHEMA_NAME = "public";
    private static final String INIT_TABLE_NAME = "mock_data";

    protected PostgresDataProviderTest(Provider provider) {
        super(provider);
    }

    @DisplayName("Test of getSchemas() method with correct DataConnection")
    @ParameterizedTest(name = "CORRECT ARGUMENTS")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#getSchemas_ok")
    @Sql(path = "scripts/postgres/postgres_basic_types.sql",
            restorePath = "scripts/postgres/drop.sql")
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
    @Sql(path = "scripts/postgres/postgres_basic_types.sql",
            restorePath = "scripts/postgres/drop.sql")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#getSchema_ok")
    public void getSchema_ok(DataFrame expected) {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getSchema(connection,
                INIT_SCHEMA_NAME, INIT_TABLE_NAME));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @Disabled
    @Test
    public void getSchema_notOk() {
        // method probably should throw something when bad input
    }

    @DisplayName("Output support for postgresql array type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_arrayType_ok")
    @Sql(path = "scripts/postgres/postgres_array.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_arrayType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for postgresql basic types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_basicTypes_ok")
    @Sql(path = "scripts/postgres/postgres_basic_types.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_basicTypes_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for postgresql bit string type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_bitStringType_ok")
    @Sql(path = "scripts/postgres/postgres_bit_string.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_bitStringType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

// bytea types - supported as string if in sql query use encode()
    @DisplayName("Output support for postgresql bytea type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_byteAType_ok")
    @Sql(path = "scripts/postgres/postgres_bytea.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_byteAType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for postgresql composite custom type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_compositeType_ok")
    @Sql(path = "scripts/postgres/postgres_composite.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_compositeType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for postgresql date, time, timestamp, interval types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_dateTypes_ok")
    @Sql(path = "scripts/postgres/postgres_dates.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_dateTypes_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for postgresql jsonb")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_jsonbType_ok")
    @Sql(path = "scripts/postgres/postgres_jsonb.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_jsonbType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for postgresql numeric, real, double precision, bigint")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_numericType_ok")
    @Sql(path = "scripts/postgres/postgres_numeric.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_numericType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for serial type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_serialType_ok")
    @Sql(path = "scripts/postgres/postgres_serial.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_serialType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for uuid type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_uuidType_ok")
    @Sql(path = "scripts/postgres/postgres_uuid.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_uuidType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for xml type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkOutputDataFrame_xmlType_ok")
    @Sql(path = "scripts/postgres/postgres_xml.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkOutputDataFrame_xmlType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource({"grok_connect.providers.arguments_provider.CommonObjectsMother#checkParameterSupport_ok",
            "grok_connect.providers.arguments_provider.CommonObjectsMother#checkMultipleParametersSupport_ok",
            "grok_connect.providers.arguments_provider.CommonObjectsMother#checkListParameterSupport_ok",
            "grok_connect.providers.arguments_provider.CommonObjectsMother#checkRegexSupport_ok"})
    @Sql(path = "scripts/postgres/postgres_basic_types.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkParameterSupport_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support for datetime")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.CommonObjectsMother#checkDatesParameterSupport_ok")
    @Sql(path = "scripts/postgres/postgres_dates_patterns.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkDatesParameterSupport_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Postgres Null safety")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.CommonObjectsMother#checkNullSupport_ok")
    @Sql(path = "scripts/postgres/postgres_null.sql",
            restorePath = "scripts/postgres/drop.sql")
    public void checkNullSupport_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall) {
        funcCall.func.connection = connection;
        Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
    }

    @DisplayName("Postgresql operators support")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.PostgresObjectsMother#checkPostgresOperatorsSupport_ok")
    @Sql(path = "scripts/postgres/postgres_operators.sql", restorePath = "scripts/postgres/drop.sql")
    public void checkPostgresOperatorsSupport_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }
}
