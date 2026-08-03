package grok_connect.connectors_info;

/**
 * RFC 8693 (OAuth 2.0 Token Exchange) descriptor. Datlas always POSTs
 * {@code application/x-www-form-urlencoded} with
 * {@code grant_type=urn:ietf:params:oauth:grant-type:token-exchange}.
 *
 * <p>Currently covers unauthenticated exchanges where the subject-token JWT
 * is alone sufficient — Databricks {@code /oidc/v1/token} is the only
 * caller today. Client credentials authentication on the exchange request
 * (Basic auth, mTLS, JWT-bearer, signed assertions) is intentionally
 * <strong>not modelled</strong>. Providers that need it (GCP STS workload
 * identity federation, Azure AD on-behalf-of when the resource requires
 * confidential-client auth, AWS STS XML dialect) will need additional
 * fields here or a sibling spec type. RFC 8693 also leaves room for
 * {@code audience}, {@code resource}, and {@code requested_token_type} —
 * add those when a provider actually needs them rather than speculatively.
 *
 * <p>{@link #endpointTemplate} supports plain {@code {name}} substitution from
 * {@code conn.parameters}. Path/query positions are URL-encoded by Datlas;
 * authority positions (host[:port]) are not.
 */
public class TokenExchangeSpec {
    /** IdP flavour this exchange applies to (e.g. {@code oidc}, {@code azure}). */
    public String flavour;

    /** Endpoint URL. {@code {paramName}} placeholders are substituted from
     *  {@code conn.parameters} (URL-encoded for path/query). */
    public String endpointTemplate;

    /** Field on the IdP response to send as {@code subject_token}.
     *  Defaults to {@code access_token} when null. */
    public String subjectTokenField;

    /** RFC 8693 token-type URN, e.g. {@code urn:ietf:params:oauth:token-type:jwt}. */
    public String subjectTokenType;

    /** Optional {@code scope} field on the exchange request. */
    public String scope;

    // RFC 8693 leaves room for `audience`, `resource`, `requested_token_type`,
    // and client-credentials authn — add when a provider actually needs them.

    public TokenExchangeSpec endpointTemplate(String t) { this.endpointTemplate = t; return this; }
    public TokenExchangeSpec subjectTokenField(String f) { this.subjectTokenField = f; return this; }
    public TokenExchangeSpec subjectTokenType(String t) { this.subjectTokenType = t; return this; }
    public TokenExchangeSpec scope(String s) { this.scope = s; return this; }
}
