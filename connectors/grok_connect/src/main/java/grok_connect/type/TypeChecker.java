package grok_connect.type;

public interface TypeChecker {
    String DEFAULT_LOG_MESSAGE = "method isSupported was called with arguments: type {}, typeName {}, "
            + "precision {}, scale {}";

    boolean isSupported(int type, String typeName, int precision, int scale);
}
