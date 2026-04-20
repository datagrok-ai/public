package grok_connect.connectors_info;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Lazy-consent OAuth/OpenID descriptor declared by a data provider.
 *
 * Scopes are keyed by IdP flavour ({@code oidc} or {@code azure}) so
 * a single connector can request different scope shapes depending on
 * which OpenID Provider the Datagrok deployment is federated with.
 * Azure AD uses the {@code <resource-uri>/.default} convention; every
 * other OIDC provider (Okta, Auth0, Keycloak, Ping, Google Workspace)
 * uses plain scopes.
 *
 * The optional {@code tokenProperty} field names the JDBC driver
 * property into which the resolved access token is written. Most
 * providers override {@code setOAuthProperty} instead, so this field
 * is purely informational — it lets Datlas or tooling display the
 * sink name without touching provider code.
 *
 * See {@code core/docs/design/GENERALIZED_OAUTH_CONNECTORS.md}.
 */
public class OAuthSpec {
    public final Map<String, List<String>> scopes = new LinkedHashMap<>();
    public String tokenProperty;

    public OAuthSpec scopes(String flavour, List<String> list) {
        scopes.put(flavour, list);
        return this;
    }

    public OAuthSpec tokenProperty(String name) {
        this.tokenProperty = name;
        return this;
    }
}
