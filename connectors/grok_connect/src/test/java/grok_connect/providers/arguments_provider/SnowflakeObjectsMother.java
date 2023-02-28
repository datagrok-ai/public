package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.DateParser;
import grok_connect.providers.utils.FuncCallBuilder;
import grok_connect.providers.utils.Parser;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.*;

import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.util.stream.Stream;

public class SnowflakeObjectsMother {
    public static Stream<Arguments> checkOutputDataFrame_dateTypes_ok() {
        Parser parser = new DateParser();
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd",
                        "2011-10-29")}), "D")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("yyyy-MM-dd HH:mm:ss.SS",
                        "2020-03-12 01:02:03.123")}), "T")
                .setColumn(new DateTimeColumn(new Double[]{parser.parseDateToDouble("HH:mm:ss.SSS X",
                        "04:05:06.789 +00:00")}), "TM")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM dates_table");
        return Stream.of(Arguments.of(Named.of("DATE TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_numericTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new BigIntColumn(new String[]{"99999999999999999999999999999999999999",
                        "-99999999999999999999999999999999999999", "9", "10000"}), "NUM")
                .setColumn(new FloatColumn(new Float[]{255.5f, 15.5f, 255.2f, 44444.4f}), "NUM10")
                .setColumn(new FloatColumn(new Float[]{Float.NaN, Float.POSITIVE_INFINITY,
                        Float.NEGATIVE_INFINITY, 1.234E+2f}), "FL1")
                .setColumn(new FloatColumn(new Float[]{0.2f, 1.34f, 317f, 124.412f}), "FL2")
                .setColumn(new IntColumn(new Integer[]{100, 51200, 7780, 12}), "I")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM test_number");
        return Stream.of(Arguments.of(Named.of("NUMERIC TYPES SUPPORT", funcCall), expected));
    }

    public static Stream<Arguments> checkOutputDataFrame_binaryType_ok() {
        String expectedHexRepr = String.format("%040x", new BigInteger(1, "Datagrok".getBytes(StandardCharsets.UTF_8)));
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{expectedHexRepr}), "B")
                .build();
        FuncCall funcCall = FuncCallBuilder.fromQuery("SELECT * FROM binary_table");
        return Stream.of(Arguments.of(Named.of("BINARY TYPE SUPPORT", funcCall), expected));
    }
}
