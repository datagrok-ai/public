package grok_connect.type;

public interface TypeChecker {
    boolean isSupported(int type, String typeName, int precision, int scale);
}
