package grok_connect.providers.arguments_provider;

import grok_connect.connectors_info.FuncCall;
import grok_connect.providers.utils.FuncCallBuilder;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.DataFrame;
import java.util.stream.Stream;

public class NeptuneObjectsMother {
    public static Stream<Arguments> someTest() {
        FuncCall funcCall = FuncCallBuilder.fromQuery("match (b:batch {source_system:'CMC'}) -[:batch_has_result_set]->" +
                "(ars:assay_result_set) return distinct b");
        return Stream.of(
                Arguments.of(Named.of("Some test", funcCall), new DataFrame()));
    }
}
