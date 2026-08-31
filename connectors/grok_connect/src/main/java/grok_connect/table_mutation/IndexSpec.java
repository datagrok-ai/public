package grok_connect.table_mutation;

import java.util.List;

/** Index definition for CreateTable; plain Gson class, no #type. */
public class IndexSpec {
    public String name; // auto-generated (ix_<table>_<col>...) when absent
    public List<String> columns;
    public boolean unique = false;
}
