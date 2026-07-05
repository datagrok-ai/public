package grok_connect.table_mutation;

import java.util.List;
import java.util.Map;

public class MutationResult {
    public int affectedRows;
    public List<PerStatementResult> perStatement;
    public List<Map<String, Object>> generatedKeys;
    public List<RowError> errors;
    // batch counters (domain schemas §5.6); null = provider cannot distinguish (filled by WO-6)
    public Integer inserted;
    public Integer updated;
    public Integer skipped;
    public Integer errorCount;
}
