package grok_connect.connectors_info;

import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.*;


public class Credentials {
    private static final String LOGIN = "login";
    private static final String PASSWORD = "password";

    public Credentials() {
    }

    public Map<String, Object> parameters = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);

    public String getLogin() { return (String)parameters.get(LOGIN); }
    public String getPassword () { return (String)parameters.get(PASSWORD); }

    public String getPasswordUrlEncoded () {
        try {
            return URLEncoder.encode(getPassword(), "UTF-8");
        } catch (UnsupportedEncodingException ex) {
            throw new RuntimeException(ex.getCause());
        }
    }
}
