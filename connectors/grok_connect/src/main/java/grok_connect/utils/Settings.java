package grok_connect.utils;

public class Settings {

    private static volatile Settings instance;

    public boolean debug;
    public int connectionPoolTimerRate;

    public Settings() {
        synchronized(Settings.class) {
            instance = this;
        }
    }

    public static Settings getInstance() {
        return instance;
    }

}
