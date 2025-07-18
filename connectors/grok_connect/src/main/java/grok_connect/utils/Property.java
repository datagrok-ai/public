package grok_connect.utils;

import java.util.List;

public class Property {
    public static final String BOOL_TYPE = "bool";
    public static final String INT_TYPE = "int";
    public static final String FLOAT_TYPE = "double";
    public static final String STRING_TYPE = "string";
    public static final String DATETIME_TYPE = "datetime";
    public static final String MAP_TYPE = "map";

    public String name;
    public String propertyType;
    public String propertySubType;
    public String description;
    public String category;
    public String format;
    public List<String> choices;
    public Prop info;

    public Property(String propertyType, String name) {
        this.name = name;
        this.propertyType = propertyType;
    }

    public Property(String propertyType, String name, Prop info) {
        this.name = name;
        this.propertyType = propertyType;
        this.info = info;
    }

    public Property(String propertyType, String name, String description) {
        this.name = name;
        this.propertyType = propertyType;
        this.description = description;
    }

    public Property(String propertyType, String name, String description, String category) {
        this.name = name;
        this.propertyType = propertyType;
        this.description = description;
        this.category = category;
    }

    public Property(String propertyType, String name, String description, Prop info) {
        this.name = name;
        this.propertyType = propertyType;
        this.description = description;
        this.info = info;
    }

    public Property(String propertyType, String name, String description, String category, Prop info) {
        this.name = name;
        this.propertyType = propertyType;
        this.description = description;
        this.category = category;
        this.info = info;
    }

    public Property(String propertyType, String name, String description, String category, Prop info, String format) {
        this.name = name;
        this.propertyType = propertyType;
        this.description = description;
        this.category = category;
        this.info = info;
        this.format = format;
    }
}
