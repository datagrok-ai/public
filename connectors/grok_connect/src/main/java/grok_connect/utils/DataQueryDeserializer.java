package grok_connect.utils;

import com.google.gson.*;
import grok_connect.connectors_info.DataQuery;
import grok_connect.table_query.TableQuery;
import java.lang.reflect.Type;

public class DataQueryDeserializer implements JsonDeserializer<DataQuery> {
    private static final Gson gson = new Gson();

    @Override
    public DataQuery deserialize(JsonElement jsonElement, Type type, JsonDeserializationContext jsonDeserializationContext) throws JsonParseException {
        String myType = jsonElement.getAsJsonObject().get("#type").getAsString();
        switch (myType) {
            case "TableQuery": return jsonDeserializationContext.deserialize(jsonElement, TableQuery.class);
            case "DataQuery": return gson.fromJson(jsonElement, DataQuery.class);
            default: return null;
        }
    }
}
