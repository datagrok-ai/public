package grok_connect.providers.arguments_provider;

import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.params.provider.Arguments;
import serialization.BigIntColumn;
import serialization.DataFrame;
import serialization.IntColumn;
import serialization.StringColumn;

import java.util.stream.Stream;

public class VerticaObjectsMother {
    private static final Parser parser = new DateParser();

    public static Stream<Arguments> getSchemas_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(8)
                .setColumn(new StringColumn(new String[] {"v_internal", "v_func", "store",
                                "v_catalog", "online_sales", "v_txtindex", "v_monitor", "public"}),
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
        String schema = "public";
        String table = "MOCK_DATA";
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
                .setColumn(new StringColumn(), fourthColumnName, new String[] {"int", "varchar(50)",
                        "varchar(50)", "varchar(50)", "varchar(50)", "varchar(50)",
                        "boolean", "varchar(50)", "date", "numeric(5,2)"})
                .setColumn(new BigIntColumn(), fifthColumnName, new String[] {"0", "0", "0", "0", "0",
                        "0", "0", "0", "0", "0"})
                .build();
        return Stream.of(Arguments.of(expected));
    }
}
