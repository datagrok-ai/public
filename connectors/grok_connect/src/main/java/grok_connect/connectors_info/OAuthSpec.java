package grok_connect.connectors_info;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Lazy-consent OAuth/OpenID descriptor declared by a data provider.
 *
 * <p>Scopes are keyed by IdP flavour ({@code oidc} or {@code azure}) so a
 * single connector can request different scope shapes depending on which
 * OpenID Provider the Datagrok deployment is federated with. Azure AD uses
 * the {@code <resource-uri>/.default} convention; every other OIDC provider
 * (Okta, Auth0, Keycloak, Ping, Google Workspace) uses plain scopes.
 *
 * <p>The flavour for a given connection is resolved on Datlas: an explicit
 * {@code oauthFlavour} parameter on the connection wins; otherwise the first
 * matching {@link HostSuffixRule} on {@code flavourParam} applies; otherwise
 * {@link #flavourDefault} (default {@code "oidc"}).
 *
 * <p>Optional {@link TokenExchangeSpec} entries describe RFC 8693
 * token-exchanges that Datlas performs against an upstream endpoint after
 * the IdP token response (e.g. swap an Azure AD {@code id_token} for a
 * Databricks workspace token via {@code /oidc/v1/token}). Each entry is
 * keyed by its own {@code flavour} field.
 */
public class OAuthSpec {
    public final Map<String, List<String>> scopes = new LinkedHashMap<>();
    public final List<TokenExchangeSpec> tokenExchange = new ArrayList<>();
    public String flavourParam;
    public final List<HostSuffixRule> flavourRules = new ArrayList<>();
    public String flavourDefault = "oidc";

    public OAuthSpec scopes(String flavour, List<String> list) {
        scopes.put(flavour, list);
        return this;
    }

    public OAuthSpec tokenExchange(String flavour, TokenExchangeSpec spec) {
        spec.flavour = flavour;
        tokenExchange.add(spec);
        return this;
    }

    /**
     * Reads the host from {@code conn.parameters[paramName]} and matches the
     * (flavour, suffix) pairs in order; first wins. Suffixes are matched
     * case-insensitively against the host portion of the parameter value.
     */
    public OAuthSpec flavourFromHostSuffix(String paramName, String... flavourSuffixPairs) {
        if (flavourSuffixPairs.length % 2 != 0)
            throw new IllegalArgumentException(
                    "flavourFromHostSuffix requires (flavour, suffix) pairs");
        this.flavourParam = paramName;
        for (int i = 0; i < flavourSuffixPairs.length; i += 2) {
            HostSuffixRule rule = new HostSuffixRule();
            rule.flavour = flavourSuffixPairs[i];
            rule.suffix = flavourSuffixPairs[i + 1];
            flavourRules.add(rule);
        }
        return this;
    }

    public OAuthSpec flavourDefault(String flavour) {
        this.flavourDefault = flavour;
        return this;
    }
}
