package grok_connect.providers.arguments_provider;

import grok_connect.providers.utils.DataFrameBuilder;
import grok_connect.providers.utils.FuncCallBuilder;
import org.junit.jupiter.api.Named;
import org.junit.jupiter.params.provider.Arguments;
import serialization.BoolColumn;
import serialization.DataFrame;
import serialization.FloatColumn;
import serialization.IntColumn;
import serialization.StringColumn;
import java.util.stream.Stream;

public class MongoDbObjectsMother {
    public static Stream<Arguments> checkOutputAllTypes_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(15)
                .setColumn(new StringColumn(new String[]{"642d3cf9675463d439e9d23d", null, null,
                        "642d3cf9c355b25fe3f3268a", null, null, "642d3cf936ac9251742be958", null, null,
                        "642d3cf975fe1b1ac7e7933a", null, null, "642d3cf9abd475f558e55197", null, null,}),
                        "map._id")
                .setColumn(new IntColumn(new Integer[]{0, -2147483648, -2147483648, 1, -2147483648, -2147483648,
                        2, -2147483648, -2147483648, 3, -2147483648, -2147483648, 4, -2147483648, -2147483648}),
                        "map.index")
                .setColumn(new StringColumn(new String[]{"1606b22d-f927-418c-a5aa-f03f14f9ce48", null, null,
                                "9a45da61-6f43-4512-9d2e-9dca26c4ff56", null, null,
                                "b92ad6a7-4d2f-4646-a42d-538085e29935", null, null,
                                "218bac3c-ddff-44fe-8429-192dfb235982", null, null,
                                "a8b57021-1203-4a1c-b907-1c2fcd06519e", null, null}),
                        "map.guid")
                .setColumn(new BoolColumn(new Boolean[]{false, false, false, true, false, false,
                        false, false, false, false, false,
                        false, true, false, false}), "map.isActive")
                .setColumn(new StringColumn(new String[]{"$2,679.87", null, null, "$3,609.65", null, null, "$2,155.94",
                        null, null,"$1,108.98", null, null,
                        "$3,429.22", null, null,}), "map.balance")
                .setColumn(new StringColumn(new String[]{"http://placehold.it/32x32", null, null,
                        "http://placehold.it/32x32", null, null, "http://placehold.it/32x32", null, null,
                        "http://placehold.it/32x32", null, null,
                        "http://placehold.it/32x32", null, null,}), "map.picture")
                .setColumn(new IntColumn(new Integer[]{28, -2147483648, -2147483648, 25, -2147483648, -2147483648, 31,
                        -2147483648, -2147483648, 25, -2147483648, -2147483648, 33, -2147483648, -2147483648,}),
                        "map.age")
                .setColumn(new FloatColumn(new Float[]{47.66397f, 2.6789344063684636e-34f, 2.6789344063684636e-34f,
                        82.178986f, 2.6789344063684636e-34f, 2.6789344063684636e-34f,
                        -65.5806f, 2.6789344063684636e-34f, 2.6789344063684636e-34f,
                        -2.464041f, 2.6789344063684636e-34f, 2.6789344063684636e-34f, 82.45335f,
                        2.6789344063684636e-34f, 2.6789344063684636e-34f, }), "map.latitude")
                .setColumn(new FloatColumn(new Float[]{-86.55363f, 2.6789344063684636e-34f, 2.6789344063684636e-34f,
                        -117.3157f, 2.6789344063684636e-34f, 2.6789344063684636e-34f,
                        -31.680079f, 2.6789344063684636e-34f, 2.6789344063684636e-34f,
                        -76.23148f, 2.6789344063684636e-34f, 2.6789344063684636e-34f,
                        74.44026f, 2.6789344063684636e-34f, 2.6789344063684636e-34f, }), "map.longitude")
                .setColumn(new StringColumn(new String[]{"[esse, dolor, ex, excepteur, consectetur, et, ad]", null, null,
                "[mollit, laborum, laborum, do, Lorem, dolore, minim]", null, null,
                                "[irure, ad, irure, occaecat, esse, amet, excepteur]", null, null,
                "[do, anim, labore, dolore, consequat, occaecat, ad]", null, null,
                                "[deserunt, dolore, ad, eu, minim, deserunt, sint]", null, null,}),
                        "map.tags")
                .setColumn(new IntColumn(new Integer[] {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2}),
                        "map.friends.id")
                .setColumn(new StringColumn(new String[]{"Mccoy Lindsey", "Joanna Lynch",
                        "Evangeline Combs", "Anita Nielsen", "Campbell Hatfield", "Mccormick Simon",
                        "Denise Curry", "Joyner Soto", "Chapman Baldwin", "Jerri Poole", "Minnie Hood",
                        "Talley Moody", "Wilkins Vaughan", "Bettie Carson", "House Pope"}), "map.friends.name")
                .build();
        return Stream.of(Arguments.of(Named.of("ALL TYPES",
                FuncCallBuilder.fromQuery("db.mocks.find()")), expected));
    }

    public static Stream<Arguments> checkStringReturn_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(4)
                .setColumn(new StringColumn(new String[]{"admin", "config", "local", "test"}),
                        "DATABASE_NAME")
                .build();
        return Stream.of(Arguments.of(Named.of("STRING RETURN TYPE",
                FuncCallBuilder.fromQuery("show dbs")), expected));
    }

    public static Stream<Arguments> checkOneDocument_ok() {
        DataFrame expected = DataFrameBuilder.getBuilder()
                .setRowCount(1)
                .setColumn(new StringColumn(new String[]{"642d3cf9675463d439e9d44d"}), "map._id")
                .setColumn(new StringColumn(new String[]{"Gomez"}),
                        "map.name")
                .setColumn(new BoolColumn(new Boolean[]{true}), "map.isActive")
                .setColumn(new StringColumn(new String[]{"$1,316.05"}),
                        "map.balance")
                .build();
        return Stream.of(Arguments.of(Named.of("DOCUMENT RETURN TYPE",
                FuncCallBuilder.fromQuery("db.one_line.find()")), expected));
    }

    public static Stream<Arguments> checkNullResult_ok() {
        return Stream.of(Arguments.of(Named.of("NULL RETURN",
                FuncCallBuilder.fromQuery("db.some_foo_bar.find()")), new DataFrame()));
    }
}
