package grok_connect.table_mutation;

public class DropIndex extends DdlMutation {
    public String indexName;

    public DropIndex() {
        type = "DropIndex";
    }
}
