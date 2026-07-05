package grok_connect.table_mutation;

import java.sql.BatchUpdateException;
import java.sql.SQLException;

import grok_connect.providers.JdbcDataProvider;

/**
 * Maps a JDBC {@link SQLException}'s SQLState to the structured per-row error {@code code}
 * used by the domain-schemas batch contract (connector-writes WO-6): {@code unique}, {@code fk},
 * {@code notnull}, {@code check}, {@code value}, {@code conflict}, else {@code db-error}. The
 * provider supplies driver-specific column/constraint extraction via {@link JdbcDataProvider}
 * hooks (Postgres reads {@code getServerErrorMessage()}); other drivers get code + message only.
 */
public class SqlStateMapper {
    public static final String DB_ERROR = "db-error";

    private SqlStateMapper() {}

    /** SQLState prefix map; {@code null}/unknown state falls back to {@link #DB_ERROR}. */
    public static String code(String sqlState) {
        if (sqlState == null || sqlState.length() < 2)
            return DB_ERROR;
        switch (sqlState) {
            case "23505": return "unique";
            case "23503": return "fk";
            case "23502": return "notnull";
            case "23514": return "check";
            case "40001":
            case "40P01": return "conflict";
            default: break;
        }
        if (sqlState.startsWith("22"))
            return "value";
        return DB_ERROR;
    }

    /** Builds a {@link RowError} for row {@code index}, unwrapping a batch exception to the real cause. */
    public static RowError toRowError(JdbcDataProvider provider, int index, SQLException e) {
        SQLException cause = unwrap(e);
        RowError error = new RowError();
        error.index = index;
        error.code = code(cause.getSQLState());
        error.column = provider.mutationErrorColumn(cause);
        error.message = provider.mutationErrorMessage(cause);
        return error;
    }

    /** A {@link BatchUpdateException} carries the failing statement's exception as its next-exception. */
    static SQLException unwrap(SQLException e) {
        if (e instanceof BatchUpdateException && e.getNextException() != null)
            return e.getNextException();
        return e;
    }
}
