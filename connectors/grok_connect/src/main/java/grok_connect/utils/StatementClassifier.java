package grok_connect.utils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;

/**
 * Conservative first-keyword classification of raw SQL statements for the read-only query policy
 * (connector-writes WO-B13, ARCHITECTURE §6.2). A statement classifies as a READ only when its first
 * keyword — after skipping whitespace, leading semicolons and parentheses, line comments, and block
 * comments — is on the read whitelist. Everything else (including an empty or unparseable statement)
 * classifies as a write. Deliberately a read-whitelist, never a write-denylist: unknown statements are
 * refused, not allowed. Best-effort by design — a CTE such as {@code WITH x AS (INSERT ...) SELECT ...}
 * passes as a read; the driver read-only session ({@code Connection.setReadOnly}) is the backstop for
 * those where enforcement is real ({@code DataSource.readOnlySessionEnforced}).
 */
public class StatementClassifier {
    private static final Set<String> READ_KEYWORDS = new HashSet<>(Arrays.asList(
            "SELECT", "WITH", "SHOW", "EXPLAIN", "DESCRIBE", "VALUES"));

    private StatementClassifier() {
    }

    public static boolean isRead(String statement, String commentStart) {
        return READ_KEYWORDS.contains(firstKeyword(statement, commentStart));
    }

    /**
     * The first keyword of {@code statement} (upper-cased), or {@code null} when none is found —
     * {@code commentStart} is the provider's line-comment marker ({@code DataSource.commentStart});
     * block comments are always skipped.
     */
    public static String firstKeyword(String statement, String commentStart) {
        if (statement == null)
            return null;
        int i = 0;
        int n = statement.length();
        while (i < n) {
            char c = statement.charAt(i);
            if (Character.isWhitespace(c) || c == '(' || c == ';') {
                i++;
                continue;
            }
            if (GrokConnectUtil.isNotEmpty(commentStart) && statement.startsWith(commentStart, i)) {
                int eol = statement.indexOf('\n', i);
                if (eol < 0)
                    return null;
                i = eol + 1;
                continue;
            }
            if (statement.startsWith("/*", i)) {
                int end = statement.indexOf("*/", i + 2);
                if (end < 0)
                    return null;
                i = end + 2;
                continue;
            }
            break;
        }
        int start = i;
        while (i < n && Character.isLetter(statement.charAt(i)))
            i++;
        return i == start ? null : statement.substring(start, i).toUpperCase(Locale.ROOT);
    }
}
