package grok_connect.providers.arguments_provider;

import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.params.provider.Arguments;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;

import java.util.stream.Stream;

public class ImpalaObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(6)
                .setColumn(new StringColumn(new String[] {"default"}),
                        "TABLE_SCHEMA")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String schema = "default";
        String table = "MOCK_DATA";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "dat", "some_number"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"BIGINT", "STRING",
                        "STRING", "STRING", "STRING", "STRING",
                        "BOOLEAN", "STRING", "DATE", "DECIMAL(5,2)"})
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table})
                .build();
        return Stream.of(Arguments.of(expected));
    }
}
