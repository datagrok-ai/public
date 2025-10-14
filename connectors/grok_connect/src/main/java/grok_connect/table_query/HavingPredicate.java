package grok_connect.table_query;

import grok_connect.utils.GrokConnectUtil;

public class HavingPredicate extends FieldPredicate {
    public String aggType;

    public HavingPredicate(String aggType, String field, String pattern, String dataType) {
        super(field, pattern, dataType);
        this.aggType = aggType;
    }

    public String getParamName() {
        String name = super.getParamName();
        return GrokConnectUtil.isEmpty(aggType) ? name : aggType.toLowerCase() + GrokConnectUtil.capitalize(name);
    }
}
