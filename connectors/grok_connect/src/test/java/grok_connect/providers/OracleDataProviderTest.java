package grok_connect.providers;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.Providers;
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

    @DisplayName("Parameters support")
    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.data_providers.CommonObjectsMother#checkParameterSupport_ok")
    public void checkParameterSupport_ok(FuncCall funcCall, DataFrame expected) {
        funcCall.func.connection = connection;
        expected.columns.removeIf(column -> column.name.equals("bool"));
        expected.columns.forEach(column -> {
            if (column.name.equals("date")) {
                column.name = "DAT";
            } else {
                column.name = column.name.toUpperCase();
            }
        });
        DataFrame actual = Assertions.assertDoesNotThrow(() -> provider.execute(funcCall));
        Assertions.assertTrue(dataFrameComparator.isDataFramesEqual(expected, actual));
    }
}