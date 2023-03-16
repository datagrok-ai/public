package grok_connect.utils;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;
import java.util.*;

public class QueryMonitor {
    private List<String> statementIdsToCancel;
    private Multimap<String, Statement> runningStatements;
    private Multimap<String, ResultSet> runningResultSets;
    private List<String> cancelledStatementIds;

    public QueryMonitor() {
        statementIdsToCancel = Collections.synchronizedList(new ArrayList<>());
        runningStatements = ArrayListMultimap.create();
        runningResultSets = ArrayListMultimap.create();
        cancelledStatementIds = Collections.synchronizedList(new ArrayList<>());
    }

    public boolean addNewStatement(String id, Statement statement) {
        synchronized(QueryMonitor.class) {
            if (id == null)
                return true;
            if (statementIdsToCancel.contains(id)) {
                statementIdsToCancel.remove(id);
                return false;
            }
            runningStatements.put(id, statement);
            return true;
        }
    }

    public boolean addNewResultSet(String id, ResultSet rs) {
        synchronized(QueryMonitor.class) {
            if (id == null)
                return true;
            
            runningResultSets.put(id, rs);
            return true;
        }
    }

    public void removeStatement(String id) {
        if (id != null)
            runningStatements.removeAll(id);
    }

    public void cancelStatement(String id) {
        synchronized(QueryMonitor.class) {
            if (runningStatements.containsKey(id)) {
                runningStatements.get(id).forEach(s -> {
                    try {
                        s.cancel();
                        runningStatements.removeAll(id);
                        cancelledStatementIds.add(id);
                    }
                    catch (SQLException throwables) {
                        throwables.printStackTrace();
                    }
                });
            }
        }
    }

    public void cancelResultSet(String id) {
        synchronized(QueryMonitor.class) {
            if (runningResultSets.containsKey(id)) {
                runningResultSets.get(id).forEach(rs -> {
                    try {
                        if (!rs.isClosed())
                            rs.close();
                        runningResultSets.removeAll(id);
                    }
                    catch (SQLException throwables) {
                        throwables.printStackTrace();
                    }
                    
                });
            }
        }
    }

    public boolean removeResultSet(String id) {
        if (runningResultSets.containsKey(id)) {
            runningResultSets.removeAll(id);
            return true;
        } 
        return false;
    }

    public boolean checkCancelledId(String id) {
        return cancelledStatementIds.remove(id);
    }
}
