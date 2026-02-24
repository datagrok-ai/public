package grok_connect.providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.NamedArgumentConverter;
import grok_connect.providers.utils.Provider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.converter.ConvertWith;
import org.junit.jupiter.params.provider.MethodSource;
import serialization.DataFrame;

class OracleDataProviderTest extends ContainerizedProviderBaseTest {
    protected OracleDataProviderTest(Provider provider) {
        super(provider);
    }

    @DisplayName("Test of getSchemas() method with correct DataConnection")
    @ParameterizedTest(name = "CORRECT ARGUMENTS")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#getSchemas_ok")
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
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#getSchema_ok")
    public void getSchema_ok(DataFrame expected) {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getSchema(connection,
                "DATAGROK", "MOCK_DATA", false));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource({"grok_connect.providers.arguments_provider.OracleObjectsMother#checkParameterSupport_ok",
            "grok_connect.providers.arguments_provider.OracleObjectsMother#checkMultipleParametersSupport_ok",
            "grok_connect.providers.arguments_provider.OracleObjectsMother#checkListParameterSupport_ok",
            "grok_connect.providers.arguments_provider.OracleObjectsMother#checkRegexSupport_ok"})
    public void checkParameterSupport_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        prepareDataFrame(expected);
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support for datetime")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkDatesParameterSupport_ok")
    public void checkDatesParameterSupport_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        expected.getColumns().forEach(column -> column.setName(column.getName().toUpperCase()));
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for xml type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_xmlType_ok")
    public void checkOutputDataFrame_xmlType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for character types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_characterTypes_ok")
    public void checkOutputDataFrame_characterTypes_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for date types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_dateTypes_ok")
    public void checkOutputDataFrame_dateTypes_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for json type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_jsonType_ok")
    public void checkOutputDataFrame_jsonType_ok(@ConvertWith(NamedArgumentConverter.class) FuncCall funcCall, DataFrame expected) {
        // resultSet.getObject(1, OracleJsonObject.class)
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for uri type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_uriType_ok")
    public void checkOutputDataFrame_uriType_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        // should be used only with .getURI()
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for varray type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_varrayType_ok")
    public void checkOutputDataFrame_varrayType_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for numeric types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_numericTypes_ok")
    public void checkOutputDataFrame_numericTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for lobs types")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.OracleObjectsMother#checkOutputDataFrame_lobsTypes_ok")
    public void checkOutputDataFrame_lobsTypes_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Oracle null safety")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.CommonObjectsMother#checkNullSupport_ok")
    public void checkNullSupport_ok(@ConvertWith(NamedArgumentConverter.class)FuncCall funcCall) {
        funcCall.func.connection = connection;
        Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
    }

    private void prepareDataFrame(DataFrame dataFrame) {
        // in order to save time reuse some common's
        dataFrame.removeColumn("bool"); // oracle doesn't have boolean type
        dataFrame.getColumns().forEach(column -> { // all columns name stored in uppercase
            if (column.getName().equals("date")) { // 'date' is reserved word in oracle, so use 'dat' for column name
                column.setName("DAT");
            } else {
                column.setName(column.getName().toUpperCase());
            }
        });
    }
}
