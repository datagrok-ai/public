package grok_connect.utils;

public class Settings {
    public boolean debug;
    public int connectionPoolTimerRate;
    public int connectionPoolMaximumPoolSize;
    public int connectionPoolIdleTimeout;

    public Settings(boolean debug, int connectionPoolTimerRate, int connectionPoolMaximumPoolSize, int connectionPoolIdleTimeout) {
        this.debug = debug;
        this.connectionPoolTimerRate = connectionPoolTimerRate;
        this.connectionPoolMaximumPoolSize = connectionPoolMaximumPoolSize;
        this.connectionPoolIdleTimeout = connectionPoolIdleTimeout;
    }
}
