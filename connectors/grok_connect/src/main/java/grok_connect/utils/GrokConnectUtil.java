package grok_connect.utils;

public class GrokConnectUtil {
    public static boolean isNotEmpty(String s) {
        return s != null && !s.isEmpty();
    }

    public static boolean isEmpty(String s) {
        return s == null || s.isEmpty();
    }
}
