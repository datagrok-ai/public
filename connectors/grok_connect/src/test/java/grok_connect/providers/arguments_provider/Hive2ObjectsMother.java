package grok_connect.providers.arguments_provider;

import grok_connect.providers.utils.DataFrameBuilder;
import org.junit.jupiter.params.provider.Arguments;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;

import java.util.stream.Stream;

public class Hive2ObjectsMother {
    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(2)
                .setColumn(new StringColumn(new String[] {"datagrok", "default"}), "table_schema")
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
                .setColumn(new StringColumn(), thirdColumnName, new String[] {"bool", "country", "date",
                        "email", "first_name", "gender", "id",
                         "ip_address",  "last_name", "some_number"})
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"boolean", "varchar(50)",
                        "date", "varchar(50)", "varchar(50)", "varchar(50)",
                        "bigint", "varchar(50)", "varchar(50)", "float"})
                .setColumn(new IntColumn(), fifthColumnName, new Integer[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
                .build();
        return Stream.of(Arguments.of(expected));
    }
}
