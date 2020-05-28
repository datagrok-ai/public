package grok_connect.utils;

import com.google.gson.*;
import java.lang.reflect.*;
import java.util.List;


public class PropertyAdapter implements JsonSerializer<Property> {

    private final Gson gson = new GsonBuilder().create();

    @Override
    public JsonElement serialize(Property src, Type typeOfSrc, JsonSerializationContext context) {
        JsonObject obj = new JsonObject();

        obj.addProperty("#type", "Property");
        obj.addProperty("name", src.name);
        obj.addProperty("propertyType", src.propertyType + (src.propertySubType == null ? "" : "<" + src.propertySubType + ">"));
        obj.addProperty("propertySubType", src.propertySubType);
        obj.addProperty("description", src.description);

        obj.add("choices", src.choices != null ? gson.toJsonTree(src.choices, List.class) : null);
        obj.add("info", src.info != null ? gson.toJsonTree(src.info, Prop.class) : null);

        return obj;
    }
}
