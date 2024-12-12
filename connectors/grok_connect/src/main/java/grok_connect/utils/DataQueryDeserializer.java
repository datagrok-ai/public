package grok_connect.utils;

import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonParseException;
import grok_connect.connectors_info.DataQuery;
import grok_connect.table_query.TableQuery;
import java.lang.reflect.Type;

public class DataQueryDeserializer implements JsonDeserializer<DataQuery> {
    @Override
    public DataQuery deserialize(JsonElement jsonElement, Type type, JsonDeserializationContext jsonDeserializationContext) throws JsonParseException {
        String myType = jsonElement.getAsJsonObject().get("#type").getAsString();
        switch (myType) {
            case "TableQuery": return jsonDeserializationContext.deserialize(jsonElement, TableQuery.class);
            case "DataQuery": return jsonDeserializationContext.deserialize(jsonElement, DataQuery.class);
            default: return null;
        }
    }
}
