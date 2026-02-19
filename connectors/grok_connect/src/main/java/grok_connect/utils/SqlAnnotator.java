package grok_connect.utils;

import grok_connect.GrokConnect;
import grok_connect.connectors_info.DataQuery;
import grok_connect.providers.JdbcDataProvider;
import net.sf.jsqlparser.expression.Expression;
import net.sf.jsqlparser.expression.ExpressionVisitorAdapter;
import net.sf.jsqlparser.parser.CCJSqlParserUtil;
import net.sf.jsqlparser.schema.Column;
import net.sf.jsqlparser.schema.Table;
import net.sf.jsqlparser.statement.Statement;
import net.sf.jsqlparser.statement.Statements;
import net.sf.jsqlparser.statement.select.*;
import serialization.DataFrame;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class SqlAnnotator {

    public static void annotate(DataQuery query, DataFrame table) {
        if (table.getColumnCount() == 0)
            return;

        Statements statements;
        JdbcDataProvider provider = GrokConnect.providerManager.getByName(query.connection.dataSource);

        boolean withSquareBrackets =
                provider != null && provider.descriptor.nameBrackets.contains("[");

        try {
            statements = CCJSqlParserUtil.parseStatements(
                    query.query,
                    parser -> parser
                            .withSquareBracketQuotation(withSquareBrackets)
                            .withAllowComplexParsing(true)
                            .withTimeOut(200)
            );
        } catch (Exception e) {
            return;
        }

        if (statements.isEmpty())
            return;

        Statement stmt = statements.get(statements.size() - 1);

        PlainSelect ps = ((Select) stmt).getPlainSelect();

        Map<String, TableInfo> ctes = collectCtes(ps);

        Map<String, TableInfo> tableAliases = new LinkedHashMap<>();

        processTables(ps, tableAliases, ctes);

        List<SelectProjection> projections = new ArrayList<>();

        processSelects(ps, projections);

        Iterator<serialization.Column<?>> dfCols = table.getColumns().iterator();

        for (SelectProjection p : projections) {

            if (p.kind == ProjectionKind.STAR) {
                while (dfCols.hasNext())
                    annotateStarColumn(dfCols.next(), tableAliases);
                break;
            }

            if (p.kind == ProjectionKind.TABLE_STAR) {
                TableInfo ti = tableAliases.get(p.tableAlias);
                if (ti == null)
                    continue;

                while (dfCols.hasNext())
                    annotateColumn(ti, dfCols.next());
                continue;
            }

            if (!dfCols.hasNext())
                break;

            serialization.Column<?> dfCol = dfCols.next();

            if (p.sourceColumns.size() == 1) {
                SourceColumn sourceColumn = getSourceColumn(p);

                TableInfo ti = sourceColumn.tableAlias != null
                        ? tableAliases.get(sourceColumn.tableAlias)
                        : null;

                if (ti == null && tableAliases.size() == 1)
                    ti = tableAliases.values().iterator().next();

                if (ti != null)
                    annotateColumn(ti, dfCol, sourceColumn.name);
            }
            else {
                dfCol.setTag(Tags.DbExpression, p.expression);
                dfCol.setTag(Tags.IsDerived, "true");
            }
        }

        if (GrokConnectUtil.isNotEmpty(query.connection.getDb())) {
            String db = query.connection.getDb();
            for (serialization.Column<?> col : table.getColumns())
                if (col.getTag(Tags.Db) == null)
                    col.setTag(Tags.Db, db);
        }

        table.setTag(Tags.DataConnectionId, query.connection.id);
    }

    public static void processTables(PlainSelect ps, Map<String, TableInfo> tableAliases, Map<String, TableInfo> ctes) {
        collectTable(ps.getFromItem(), tableAliases, ctes);
        if (ps.getJoins() != null) {
            for (Join j : ps.getJoins())
                collectTable(j.getRightItem(), tableAliases, ctes);
        }
    }

    public static void processSelects(PlainSelect ps, List<SelectProjection> projections) {
        for (SelectItem<?> item : ps.getSelectItems()) {

            Expression expr = item.getExpression();

            if (expr instanceof AllTableColumns) {
                AllTableColumns atc = (AllTableColumns) expr;
                projections.add(
                        SelectProjection.tableStar(
                                normalize(atc.getTable().getName())
                        )
                );
                continue;
            }


            if (expr instanceof AllColumns) {
                projections.add(SelectProjection.star());
                continue;
            }

            String alias;

            if (item.getAlias() != null)
                alias = normalize(item.getAlias().getName());
            else if (expr instanceof Column)
                alias = normalize(((Column) expr).getColumnName());
            else
                alias = normalize(expr.toString());

            Set<String> sourceCols = new LinkedHashSet<>();
            expr.accept(new ExpressionVisitorAdapter() {
                @Override
                public void visit(Column c) {
                    sourceCols.add(c.getFullyQualifiedName());
                }
            });

            projections.add(
                    SelectProjection.explicit(
                            alias,
                            expr.toString(),
                            sourceCols
                    )
            );
        }
    }

    private static ColumnLineage resolveColumnLineageThroughSource(TableInfo ti, String colName) {
        if (ti == null)
            return null;

        if (!(ti instanceof DerivedTableInfo))
            return new ColumnLineage(ti.database, ti.table, ti.schema, colName);

        DerivedTableInfo dti = (DerivedTableInfo) ti;

        return dti.columnLineage.get(colName);
    }

    private static SourceColumn getSourceColumn(SelectProjection p) {
        String fq = p.sourceColumns.iterator().next();
        String[] parts = fq.split("\\.");

        String colName = normalize(parts[parts.length - 1]);
        String tableAlias = parts.length > 1 ? normalize(parts[parts.length - 2]) : null;
        return new SourceColumn(colName, tableAlias);
    }

    private static void collectTable(FromItem item, Map<String, TableInfo> map, Map<String, TableInfo> ctes) {
        if (item instanceof ParenthesedSelect)
            processParenthesedSelect((ParenthesedSelect) item, item.getAlias().getName(), map, ctes);
        else if (item instanceof Table) {
            Table t = (Table) item;

            String alias = normalize(t.getAlias() != null ? t.getAlias().getName() : t.getName());
            String name = normalize(t.getName());

            // If table name matches a CTE name, alias points to derived table
            if (ctes != null && ctes.containsKey(name)) {
                DerivedTableInfo cte = (DerivedTableInfo) ctes.get(name);

                // Re-alias it for this FROM alias
                map.put(alias, new DerivedTableInfo(alias, cte.columnLineage));
                return;
            }

            map.put(
                    alias,
                    new TableInfo(
                            normalize(t.getDatabase().getDatabaseName()),
                            normalize(t.getSchemaName()),
                            normalize(t.getName()),
                            alias
                    )
            );
        }
    }

    private static String normalize(String s) {
        if (s == null)
            return null;
        return s.replace("\"", "")
                .replace("[", "")
                .replace("]", "")
                .replace("`", "");
    }

    private static void annotateStarColumn(
            serialization.Column<?> col,
            Map<String, TableInfo> tables
    ) {
        if (tables.size() == 1) {
            TableInfo ti = tables.values().iterator().next();
            annotateColumn(ti, col);
        }
    }

    private static void annotateColumn(TableInfo ti, serialization.Column<?> col) {
        annotateColumn(ti, col, col.getName());
    }

    private static void annotateColumn(TableInfo ti, serialization.Column<?> col, String colName) {
        if (ti instanceof DerivedTableInfo)
            annotateUsingDerived((DerivedTableInfo) ti, col, colName);
        else {
            col.setTag(Tags.DbTable, ti.table);
            col.setTag(Tags.DbColumn, col.getName());
            if (ti.schema != null)
                col.setTag(Tags.DbSchema, ti.schema);
            if (ti.database != null)
                col.setTag(Tags.Db, ti.database);
        }
    }

    private static void annotateUsingDerived(DerivedTableInfo ti, serialization.Column<?> col, String colName) {
        ColumnLineage cl = ti.columnLineage.get(colName);
        if (cl != null) {
            col.setTag(Tags.DbTable, cl.table);
            col.setTag(Tags.DbColumn, cl.column);
            if (cl.schema != null)
                col.setTag(Tags.DbSchema, cl.schema);
            if (cl.database != null)
                col.setTag(Tags.Db, cl.database);
        }
    }

    private static Map<String, TableInfo> collectCtes(Select selectStmt) {
        Map<String, TableInfo> ctes = new LinkedHashMap<>();

        List<WithItem> withItems = selectStmt.getWithItemsList();
        if (withItems == null)
            return ctes;

        // Build in order: later CTEs can reference earlier ones
        for (WithItem wi : withItems) {
            String name = normalize(wi.getAlias().getName());

            Select s = wi.getSelect();

            if (s == null) {
                ctes.put(name, new DerivedTableInfo(name, new LinkedHashMap<>()));
                continue;
            }

            processParenthesedSelect(s, name, ctes, ctes);
        }

        return ctes;
    }

    private static void processParenthesedSelect(Select ps, String alias, Map<String, TableInfo> map, Map<String, TableInfo> ctes) {
        PlainSelect plainSelect = ps.getPlainSelect();
        Map<String, TableInfo> innerTables = new LinkedHashMap<>();
        processTables(plainSelect, innerTables, ctes);
        List<SelectProjection> innerProjections = new ArrayList<>();
        processSelects(plainSelect, innerProjections);

        Map<String, ColumnLineage> lineages = buildLineageForSelect(innerTables, innerProjections);

        map.put(alias, new DerivedTableInfo(alias, lineages));
    }

    private static Map<String, ColumnLineage> buildLineageForSelect(Map<String, TableInfo> tableAliases, List<SelectProjection> projections) {
        Map<String, ColumnLineage> lineages = new LinkedHashMap<>();

        for (SelectProjection p : projections) {

            if (p.kind == ProjectionKind.STAR || p.kind == ProjectionKind.TABLE_STAR)
                continue;

            if (p.sourceColumns == null || p.sourceColumns.size() != 1)
                continue;

            SourceColumn sourceColumn = getSourceColumn(p);

            TableInfo ti = null;

            if (sourceColumn.tableAlias != null)
                ti = tableAliases.get(sourceColumn.tableAlias);

            if (ti == null && tableAliases.size() == 1)
                ti = tableAliases.values().iterator().next();

            String outName = normalize(p.alias != null ? p.alias : sourceColumn.name);

            ColumnLineage resolved = resolveColumnLineageThroughSource(ti, sourceColumn.name);

            if (resolved != null)
                lineages.put(outName, resolved);
        }
        return lineages;
    }


    private enum ProjectionKind {
        STAR,
        TABLE_STAR,
        EXPLICIT
    }

    private static class SelectProjection {
        ProjectionKind kind;
        String alias;
        String expression;
        String tableAlias;
        Set<String> sourceColumns;

        static SelectProjection star() {
            SelectProjection p = new SelectProjection();
            p.kind = ProjectionKind.STAR;
            return p;
        }

        static SelectProjection tableStar(String alias) {
            SelectProjection p = new SelectProjection();
            p.kind = ProjectionKind.TABLE_STAR;
            p.tableAlias = alias;
            return p;
        }

        static SelectProjection explicit(
                String alias, String expr, Set<String> src
        ) {
            SelectProjection p = new SelectProjection();
            p.kind = ProjectionKind.EXPLICIT;
            p.alias = alias;
            p.expression = expr;
            p.sourceColumns = src;
            return p;
        }
    }

    private static class TableInfo {
        String database;
        String schema;
        String table;
        String alias;

        TableInfo(String d, String s, String t, String a) {
            database = d;
            schema = s;
            table = t;
            alias = a;
        }
    }

    private static class DerivedTableInfo extends TableInfo {
        Map<String, ColumnLineage> columnLineage;

        DerivedTableInfo(String a, Map<String, ColumnLineage> lineages) {
            super(null, null, null, a);
            columnLineage = lineages;
        }
    }

    private static class SourceColumn {
        String name;
        String tableAlias;

        SourceColumn(String name, String tableAlias) {
            this.name = name;
            this.tableAlias = tableAlias;
        }
    }

    private static class ColumnLineage {
        String database;
        String table;
        String schema;
        String column;

        ColumnLineage(String database, String table, String schema, String column) {
            this.database = database;
            this.table = table;
            this.schema = schema;
            this.column = column;
        }
    }
}
