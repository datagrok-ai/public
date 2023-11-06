package grok_connect.utils;

public class GrokConnectException extends RuntimeException {
    public GrokConnectException(String message) {
        super(message);
    }

    public GrokConnectException(String message, Throwable cause) {
        super(message, cause);
    }
}
