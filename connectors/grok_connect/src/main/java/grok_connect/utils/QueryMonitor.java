package grok_connect.utils;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class QueryMonitor {
    private static QueryMonitor instance;
    private final Set<String> statementIdsToCancel;
    private final Set<String> resultSetIdsToCancel;
    private final Multimap<String, Statement> runningStatements;
    private final Set<String> cancelledStatementIds;

    private QueryMonitor() {
        statementIdsToCancel = Collections.synchronizedSet(new HashSet<>());
        resultSetIdsToCancel = Collections.synchronizedSet(new HashSet<>());
        runningStatements = ArrayListMultimap.create();
        cancelledStatementIds = Collections.synchronizedSet(new HashSet<>());
    }

    public static synchronized QueryMonitor getInstance() {
        if (instance == null) {
            instance = new QueryMonitor();
        }
        return instance;
    }

    public void addNewStatement(String id, Statement statement) {
        synchronized(QueryMonitor.class) {
            if (id == null)
                return;
            if (statementIdsToCancel.contains(id)) {
                statementIdsToCancel.remove(id);
                return;
            }
            runningStatements.put(id, statement);
        }
    }

    public void addCancelledResultSet(String id) {
        resultSetIdsToCancel.add(id);
    }

    public void removeStatement(String id) {
        if (id != null)
            runningStatements.removeAll(id);
    }

    public void cancelStatement(String id) {
        if (runningStatements.containsKey(id)) {
            synchronized(QueryMonitor.class) {
                cancelledStatementIds.add(id);

                Collection<Statement> statements = runningStatements.removeAll(id);
                for (Statement s: statements) {
                    try {
                        s.cancel();
                    } catch (SQLException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
    }

    public void removeResultSet(String id) {
        resultSetIdsToCancel.remove(id);
    }

    public boolean checkCancelledId(String id) {
        return cancelledStatementIds.remove(id);
    }

    public boolean checkCancelledIdResultSet(String id) {
        return resultSetIdsToCancel.contains(id);
    }
}
