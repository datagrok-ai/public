package grok_connect.table_mutation;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import grok_connect.providers.JdbcDataProvider;
import grok_connect.utils.QueryMonitor;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Plans and executes a single DDL operation under the ARCHITECTURE §3.4.2 contract: the plan carries
 * the exact per-provider statements plus live-data pre-check counts (read-only SELECTs on the same
 * connection, BEFORE any DDL); the apply path recomputes the plan (a client's earlier dry run is
 * never trusted) and refuses destructive actions without confirmDestructive. Transaction control
 * stays with {@link MutationRunner} — commit on success, rollback on any error, always before the
 * connection returns to the pool. The mutation dryRun contract is unrelated to QueryManager's
 * debug-mode dry run for queries.
 */
public class DdlRunner {
    private static final Logger LOGGER = LoggerFactory.getLogger(DdlRunner.class);

    /** Exact per-provider statements for a DDL op — the single emission path (the WO-B6 emitters). */
    public static List<String> statements(JdbcDataProvider provider, DdlMutation m) {
        if (m instanceof CreateTable)
            return provider.createTableSql((CreateTable) m);
        if (m instanceof AlterTable)
            return Collections.singletonList(provider.alterTableSql((AlterTable) m));
        if (m instanceof CreateIndex)
            return Collections.singletonList(provider.createIndexSql((CreateIndex) m));
        if (m instanceof DropIndex)
            return Collections.singletonList(provider.dropIndexSql((DropIndex) m));
        if (m instanceof AddKey)
            return Collections.singletonList(provider.addKeySql((AddKey) m));
        if (m instanceof DropTable)
            return Collections.singletonList(provider.dropTableSql((DropTable) m));
        if (m instanceof TruncateTable)
            return Collections.singletonList(provider.truncateTableSql((TruncateTable) m));
        throw new MutationValidationException("Unsupported DDL mutation type: " + m.type);
    }

    /**
     * Statement emission + destructive classification with live pre-check COUNTs. Hard refusals —
     * cases where the DDL would fail regardless of confirmation — throw a structured validation
     * error carrying {@code {path, code, message, count}} refusals. Pre-checks are read-only and
     * run before any DDL, so a refusal leaves the transaction clean for rollback.
     */
    public static MutationPlan plan(JdbcDataProvider provider, Connection connection, DdlMutation m) throws SQLException {
        MutationPlan plan = new MutationPlan();
        plan.statements = statements(provider, m); // validates the payload before any DB round trip
        plan.transactionalDdl = provider.descriptor.supportsTransactionalDdl;
        plan.destructive = new ArrayList<>();
        if (m instanceof DropTable || m instanceof TruncateTable) {
            long rows = count(provider, connection, m, null);
            boolean drop = m instanceof DropTable;
            plan.destructive.add(action("table", drop ? "drop-table" : "truncate-table",
                    (drop ? "Dropping" : "Truncating") + " table " + m.tableName
                            + " permanently deletes " + rows + " row(s)", rows));
        }
        else if (m instanceof AlterTable) {
            AlterTable alter = (AlterTable) m;
            switch (alter.action) { // per-action payload presence is already validated by statements()
                case "dropColumn": {
                    long values = count(provider, connection, m, provider.addBrackets(alter.columnName) + " IS NOT NULL");
                    plan.destructive.add(action("columns." + alter.columnName, "drop-column",
                            "Dropping column " + alter.columnName + " discards " + values + " non-null value(s)", values));
                    break;
                }
                case "changeType": {
                    long values = count(provider, connection, m, provider.addBrackets(alter.columnName) + " IS NOT NULL");
                    plan.destructive.add(action("columns." + alter.columnName, "type-change",
                            "Changing the type of column " + alter.columnName + " converts " + values
                                    + " non-null value(s); a failing cast surfaces as a database error", values));
                    break;
                }
                case "setNullable": {
                    if (!alter.nullable) {
                        long nulls = count(provider, connection, m, provider.addBrackets(alter.columnName) + " IS NULL");
                        if (nulls > 0)
                            refuse("columns." + alter.columnName, "not-null-violation", "Column " + alter.columnName
                                    + " has " + nulls + " NULL value(s); NOT NULL cannot be applied", nulls);
                    }
                    break;
                }
                case "addColumn": {
                    if (!alter.column.nullable && alter.column.defaultValue == null) {
                        long rows = count(provider, connection, m, null);
                        if (rows > 0)
                            refuse("columns." + alter.column.name, "required-without-default",
                                    "Cannot add required column " + alter.column.name
                                            + " without a default value: table has " + rows + " row(s)", rows);
                    }
                    break;
                }
            }
        }
        return plan;
    }

    /**
     * Recomputes the plan, refuses non-empty destructive actions without confirmDestructive, then
     * runs each statement on the caller's transaction. Returns the plan for the audit value on
     * {@link MutationResult#plan}.
     */
    public static MutationPlan execute(JdbcDataProvider provider, Connection connection, DdlMutation m, String mainCallId) throws SQLException {
        MutationPlan plan = plan(provider, connection, m);
        if (!plan.destructive.isEmpty() && !m.confirmDestructive)
            throw new MutationConfirmationRequiredException(plan, provider.descriptor.type);
        QueryMonitor queryMonitor = QueryMonitor.getInstance();
        for (String sql : plan.statements) {
            LOGGER.info("Mutation before execution: {}", sql);
            try (Statement statement = connection.createStatement()) {
                queryMonitor.addNewStatement(mainCallId, statement);
                try {
                    statement.execute(sql);
                } finally {
                    queryMonitor.removeStatement(mainCallId);
                }
            }
        }
        return plan;
    }

    private static DestructiveAction action(String path, String code, String message, long count) {
        DestructiveAction action = new DestructiveAction();
        action.path = path;
        action.code = code;
        action.message = message;
        action.count = count;
        return action;
    }

    private static void refuse(String path, String code, String message, long count) {
        throw new MutationValidationException(message, Collections.singletonList(action(path, code, message, count)));
    }

    /** Read-only live-data pre-check on the mutation connection (joins its transaction; never mutates). */
    private static long count(JdbcDataProvider provider, Connection connection, TableMutation m, String where) throws SQLException {
        String sql = "SELECT COUNT(*) FROM " + provider.mutationTableName(m) + (where == null ? "" : " WHERE " + where);
        try (Statement statement = connection.createStatement();
             ResultSet rs = statement.executeQuery(sql)) {
            rs.next();
            return rs.getLong(1);
        }
    }
}
