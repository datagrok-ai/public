package grok_connect.utils;

import com.google.gson.*;
import java.lang.reflect.*;


public class PropertyAdapter implements JsonSerializer<Property> {

    private Gson gson = new GsonBuilder().create();

    private void addNullableSubClass(JsonObject obj, Object sub, String name) {
        if (sub != null)
            obj.add(name, gson.toJsonTree(sub, sub.getClass()));
        else
            obj.add(name, null);
    }

    @Override
    public JsonElement serialize(Property src, Type typeOfSrc, JsonSerializationContext context) {
        JsonObject obj = new JsonObject();

        obj.addProperty("#type", "Property");
        obj.addProperty("name", src.name);
        obj.addProperty("propertyType",
                src.propertyType + (src.propertySubType == null ? "" : "<" + src.propertySubType + ">"));
        obj.addProperty("propertySubType", src.propertySubType);
        obj.addProperty("description", src.description);

        addNullableSubClass(obj, src.choices, "choices");
        addNullableSubClass(obj, src.info, "info");

        return obj;
    }
}
