package grok_connect.providers;

import grok_connect.connectors_info.DataQuery;
import grok_connect.connectors_info.DataSource;
import grok_connect.providers.utils.NamedArgumentConverter;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.converter.ConvertWith;
import org.junit.jupiter.params.provider.MethodSource;
import org.mockito.ArgumentMatchers;
import org.mockito.Mockito;
import org.junit.jupiter.api.Assertions;
import java.util.List;

class JdbcDataProviderTest {
    private static final String NOT_EXISTED_PARAM_NAME = "name";
    private static JdbcDataProvider jdbcDataProvider;
    private static DataQuery dataQuery;

    @BeforeAll
    public static void init() {
        jdbcDataProvider = Mockito.mock(JdbcDataProvider.class, Mockito.CALLS_REAL_METHODS);
        jdbcDataProvider.descriptor = new DataSource();
        dataQuery = Mockito.mock(DataQuery.class);
        Mockito.when(dataQuery.existsParam(ArgumentMatchers.any(String.class))).thenReturn(true);
        Mockito.when(dataQuery.existsParam(NOT_EXISTED_PARAM_NAME)).thenReturn(false);
    }

    @ParameterizedTest(name = "{index} : {0}")
    @MethodSource("grok_connect.providers.arguments_provider.JdbcObjectsMother#test_getParameterNames_ok")
    public void test_getParameterNames_ok(@ConvertWith(NamedArgumentConverter.class) String query, List<String> expectedNames,
                                          StringBuilder expectedBuffer) {
        StringBuilder queryBuffer = new StringBuilder();
        List<String> actualNames = Assertions.assertDoesNotThrow(() ->
                jdbcDataProvider.getParameterNames(query, dataQuery, queryBuffer));
        Assertions.assertEquals(expectedNames, actualNames);
        Assertions.assertEquals(expectedBuffer.toString(), queryBuffer.toString());
    }
}
