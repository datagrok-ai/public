package grok_connect.table_mutation;

import java.sql.SQLException;

/**
 * A keyed UPDATE/DELETE affected a different number of rows than its {@code expectAffected}
 * (external-bindings WRITEBACK A.2 #6). Raised inside the transaction, before commit, so the
 * runner's SQLException path rolls the whole batch back; {@link SqlStateMapper} maps
 * {@link #SQL_STATE} to the {@code affected} row-error code and carries {@link #actual}.
 */
public class AffectedRowsMismatchException extends SQLException {
    public static final String SQL_STATE = "DG001";

    public final int expected;
    public final int actual;

    public AffectedRowsMismatchException(int expected, int actual) {
        super("Expected " + expected + " affected row(s), got " + actual, SQL_STATE);
        this.expected = expected;
        this.actual = actual;
    }
}
