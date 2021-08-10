package grok_connect.utils;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class QueryMonitor {
    private List<String> statementIdsToCancel;
    private Multimap<String, Statement> runningStatements;

    public QueryMonitor() {
        statementIdsToCancel = new ArrayList<>();
        runningStatements = ArrayListMultimap.create();
    }

    public boolean addNewStatement(String id, Statement statement) {
        if (statementIdsToCancel.contains(id)) {
            statementIdsToCancel.remove(id);
            return false;
        }
        runningStatements.put(id, statement);
        return true;
    }

    public void cancelStatement(String id) {
        if (runningStatements.containsKey(id)) {
            runningStatements.get(id).forEach(s -> {
                try {
                    s.cancel();
                }
                catch (SQLException throwables) {
                    throwables.printStackTrace();
                }
            });
        }
    }
}
