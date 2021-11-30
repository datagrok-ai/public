package grok_connect.utils;

public class SettingsManager {
    private static volatile SettingsManager instance;

    public Settings settings;

    private SettingsManager() {
    }

    public void setSettings(Settings settings) {
        this.settings = settings;
    }

    public static SettingsManager getInstance() {
        SettingsManager result = instance;
        if (result != null) {
            return result;
        }
        synchronized(SettingsManager.class) {
            if (instance == null) {
                instance = new SettingsManager();
            }
            return instance;
        }
    }
}
