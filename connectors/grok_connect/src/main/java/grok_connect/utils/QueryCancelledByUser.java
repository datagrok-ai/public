package grok_connect.utils;

public class QueryCancelledByUser extends Exception {
    @Override
    public String getMessage() {
        return "QueryCancelledByUser";
    }
}
