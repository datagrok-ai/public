package grok_connect.utils;

import java.util.List;
import serialization.Column;

public class SchemeInfo {
    public List<Column> columns;
    public List<Boolean> supportedType;
    public List<Boolean> initColumn;

    public SchemeInfo(List<Column> columns, List<Boolean> supportedType, List<Boolean> initColumn){
        this.columns = columns;
        this.supportedType = supportedType;
        this.initColumn = initColumn;
    }
}
