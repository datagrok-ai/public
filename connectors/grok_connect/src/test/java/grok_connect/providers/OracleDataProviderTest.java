package grok_connect.providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.Providers;
import grok_connect.providers.utils.Sql;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import serialization.DataFrame;

class OracleDataProviderTest extends ProviderBaseTest {

    protected OracleDataProviderTest() {
        super(Providers.ORACLE);
    }

    @DisplayName("Test of getSchemas() method with correct DataConnection")
    @ParameterizedTest(name = "CORRECT ARGUMENTS")
    @MethodSource("grok_connect.providers.data_providers.OracleObjectsMother#getSchemas_ok")
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
    @MethodSource("grok_connect.providers.data_providers.OracleObjectsMother#getSchema_ok")
    public void getSchema_ok(DataFrame expected) {
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.getSchema(connection,
                "DATAGROK", "MOCK_DATA"));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource({"grok_connect.providers.data_providers.CommonObjectsMother#checkParameterSupport_ok",
            "grok_connect.providers.data_providers.OracleObjectsMother#checkMultipleParametersSupport"})
    public void checkParameterSupport_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        prepareDataFrame(expected);
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Parameters support for datetime")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.data_providers.OracleObjectsMother#checkDatesParameterSupport_ok")
    public void checkDatesParameterSupport_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        expected.columns.forEach(column -> column.name = column.name.toUpperCase());
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    @DisplayName("Output support for xml type")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.data_providers.OracleObjectsMother#checkOutputDataFrame_xmlType_ok")
    public void checkOutputDataFrame_xmlType_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }

    private void prepareDataFrame(DataFrame dataFrame) {
        // in order to save time reuse some common's
        dataFrame.columns.removeIf(column -> column.name.equals("bool")); // oracle doesn't have boolean type
        dataFrame.columns.forEach(column -> { // all columns name stored in uppercase
            if (column.name.equals("date")) { // 'date' is reserved word in oracle, so use 'dat' for column name
                column.name = "DAT";
            } else {
                column.name = column.name.toUpperCase();
            }
        });
    }
}