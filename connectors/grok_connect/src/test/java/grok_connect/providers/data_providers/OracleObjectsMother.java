package grok_connect.providers.data_providers;

import grok_connect.providers.utils.DataFrameBuilder;
import org.junit.jupiter.params.provider.Arguments;
import serialization.DataFrame;
import serialization.StringColumn;

import java.util.stream.Stream;

/**
 * Provides data for Oracle provider tests
 */
public class OracleObjectsMother {
    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"DATAGROK", "GSMADMIN_INTERNAL",
                        "LBACSYS"}), "TABLE_SCHEMA")
                .build();
        return Stream.of(Arguments.of(expected));
    }
}
