package grok_connect.utils;

import com.google.gson.annotations.*;


public class Prop {
    public String editor;
    public boolean ignore;

    public Prop() {}

    public Prop(String editor) {
        this.editor = editor;
    }

    public Prop(String editor, boolean ignore) {
        this.editor = editor;
        this.ignore = ignore;
    }
}
