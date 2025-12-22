package grok_connect.utils;

public class GrokConnectException extends Exception {
    public GrokConnectException(String message) {
        super(message);
    }

    public GrokConnectException(String message, Throwable cause) {
        super(message, cause);
    }

    public GrokConnectException(Throwable cause) {
        super(cause);
    }
}
