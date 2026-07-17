package grok_connect.table_mutation;

import java.util.List;

public class CreateIndex extends DdlMutation {
    public String indexName; // auto-generated (ix_<table>_<col>...) when absent
    public List<String> columns;
    public boolean unique;

    public CreateIndex() {
        type = "CreateIndex";
    }
}
