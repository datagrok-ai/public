package grok_connect.connectors_info;

/**
 * One ordered rule in {@link OAuthSpec#flavourRules}: if the host portion of
 * {@code conn.parameters[OAuthSpec.flavourParam]} ends with {@link #suffix}
 * (case-insensitive), the connection resolves to {@link #flavour}.
 */
public class HostSuffixRule {
    public String flavour;
    public String suffix;
}
