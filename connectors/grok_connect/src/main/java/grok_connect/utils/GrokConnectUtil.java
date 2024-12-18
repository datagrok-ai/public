package grok_connect.utils;

public class GrokConnectUtil {
    public static boolean isNotEmpty(String s) {
        return s != null && !s.isEmpty();
    }

    public static boolean isEmpty(String s) {
        return s == null || s.isEmpty();
    }

    public static String capitalize(String s) {
        return isEmpty(s) ? s : s.substring(0, 1).toUpperCase() + s.substring(1);
    }
}
