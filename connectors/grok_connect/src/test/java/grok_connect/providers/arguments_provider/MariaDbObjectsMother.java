package grok_connect.providers.arguments_provider;

import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.FuncCallBuilder;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;

import java.util.stream.Stream;

public class MariaDbObjectsMother {
    public static Stream<Arguments> checkOutputDataFrame_jsonType_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(3)
                .setColumn(new StringColumn(new String[] {"{\"key1\": \"value1\", \"key2\": \"value2\"}",
                                "{ \"phones\":[ {\"type\": \"mobile\", \"phone\": \"001001\"} , "
                                        + "{\"type\": \"fix\", \"phone\": \"002002\"} ] }", "{\"reading\": 1.230e-5}"}),
                        "json_type")
                .build();
        return Stream.of(Arguments.of(Named.of("JSON TYPE SUPPORT",
                FuncCallBuilder.fromQuery("SELECT * FROM JSON_TYPE")), expected));
    }

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[] {"datagrok", "information_schema"}),
                        "table_schema")
                .build();
        return Stream.of(Arguments.of(expected));
    }

    public static Stream<Arguments> getSchema_ok() {
        String firstColumnName = "table_schema";
        String secondColumnName = "table_name";
        String thirdColumnName = "column_name";
        String fourthColumnName = "data_type";
        String fifthColumnName = "is_view";
        String schema = "datagrok";
        String table = "mock_data";
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(10)
                .setColumn(new StringColumn(), firstColumnName, new String[] {schema, schema,
                        schema, schema, schema, schema, schema,
                        schema, schema, schema})
                .setColumn(new StringColumn(), secondColumnName, new String[] {table, table,
                        table, table, table, table,
                        table, table, table, table})
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"id", "first_name", "last_name", "email",
                        "gender", "ip_address", "bool", "country", "date", "some_number"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"bigint", "varchar",
                        "varchar", "varchar", "varchar", "varchar",
                        "tinyint", "varchar", "date", "decimal"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }
}
