package grok_connect.log;

import java.util.Map;

public class LogMessage {
    private final String level;
    private final long time;
    private final String message;
    private final String flag;
    private final Map<String, Object> params;
    private final String stackTrace;

    public LogMessage(String level, long time, String message, String flag, Map<String, Object> params, String stackTrace) {
        this.level = level;
        this.time = time;
        this.message = message;
        this.flag = flag;
        this.params = params;
        this.stackTrace = stackTrace;
    }

    public String getLevel() {
        return level;
    }

    public long getTime() {
        return time;
    }

    public String getMessage() {
        return message;
    }

    public String getFlag() {
        return flag;
    }

    public Map<String, Object> getParams() {
        return params;
    }

    public String getStackTrace() {
        return stackTrace;
    }
}
