package grok_connect.utils;

import com.google.gson.annotations.*;


public class Prop {
    public String editor;
    public boolean ignore;
    public boolean nullable = true;

    public Prop() {}

    public Prop(String editor) {
        this.editor = editor;
    }

    public Prop(String editor, boolean ignore) {
        this.editor = editor;
        this.ignore = ignore;
    }

    public Prop(String editor, boolean ignore, boolean nullable) {
        this.editor = editor;
        this.ignore = ignore;
        this.nullable = nullable;
    }
}
