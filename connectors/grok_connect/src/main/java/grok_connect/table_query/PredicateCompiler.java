package grok_connect.table_query;

import java.util.List;

import grok_connect.connectors_info.FuncParam;
import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.GrokConnectUtil;
import grok_connect.utils.PatternMatcherResult;
import serialization.Types;

/**
 * Compiles a {@link FieldPredicate} into a SQL fragment with @-named parameters through the
 * per-provider pattern converters. Shared by TableQuery.toSql and table_mutation SQL emission —
 * one WHERE-compilation path.
 */
public class PredicateCompiler {
    /** Returns the SQL fragment for the clause; appends the generated parameters to {@code out}. */
    public static String compile(FieldPredicate clause, JdbcDataProvider provider, List<FuncParam> out) {
        String paramName = clause.getParamName();
        clause.matcher.colName = provider.addBrackets(clause.field);
        if (clause instanceof HavingPredicate && GrokConnectUtil.isNotEmpty(((HavingPredicate) clause).aggType)) {
            provider.descriptor.aggregations.stream()
                    .filter((a) -> a.functionName.equals(((HavingPredicate) clause).aggType))
                    .findFirst()
                    .ifPresent(info -> clause.matcher.colName = info.dbFunctionName
                            .replaceAll("#", clause.matcher.colName));
        }
        PatternMatcherResult result;
        switch (clause.dataType) {
            case Types.NUM:
            case Types.FLOAT:
            case Types.INT:
                result = provider.numericPatternConverter(paramName, clause.dataType, clause.matcher);
                break;
            case Types.STRING:
                result = provider.stringPatternConverter(paramName, clause.matcher);
                break;
            case Types.DATE_TIME:
                result = provider.dateTimePatternConverter(paramName, clause.matcher);
                break;
            case Types.BOOL:
                result = provider.boolPatternConverter(paramName, clause.matcher);
                break;
            case Types.BIG_INT:
                result = provider.bigIntPatternConverter(paramName, clause.matcher);
                break;
            default:
                throw new UnsupportedOperationException(clause.dataType + " is not supported");
        }
        out.addAll(result.params);
        return result.query;
    }
}
