---
title: "OAuth authentication for connectors"
sidebar_label: "OAuth (BigQuery, Databricks)"
---

Datagrok supports per-user OAuth authentication for [BigQuery](bigquery.md)
and [Databricks](databricks.md). Each user signs in with their own
identity, so row-level security, dataset grants, and audit trails at the
data source apply per user instead of to a shared service account.

Consent happens lazily on first use. Datagrok stores the resulting access
and refresh tokens against the user's personal group and refreshes them
silently on every subsequent query.

## How it works

1. A user opens or queries an OAuth-enabled connection for the first time.
2. Datagrok shows a popup (on **Test** / **Run Query**) or a balloon with an
   **Authorize** button (on browse-tree expansion).
3. The user consents at the identity provider. Datagrok exchanges the code
   for access and refresh tokens and stores them encrypted.
4. Every subsequent query reuses the token, or silently refreshes it 60
   seconds before expiry. The popup reappears only after a scope change or
   a revoked grant.

## Identity provider setup

OAuth for connectors reuses the same OpenID configuration as user login
(**Settings → OpenID**). Datagrok auto-detects the flavour from the
discovery URL: `login.microsoftonline.com*` → Azure, everything else → OIDC.

Add the connector callback to your identity-provider application's redirect
URIs:

```text
https://<your-datagrok-host>/oauth/connector/callback
```

Then grant the delegated permission your connector needs:

| Provider            | Connector  | Delegated permission                                                           |
|---------------------|------------|--------------------------------------------------------------------------------|
| Azure AD            | BigQuery   | `https://bigquery.googleapis.com/.default` (requires Google [Workforce Identity Federation](https://cloud.google.com/iam/docs/workforce-identity-federation)) |
| Azure AD            | Databricks | AzureDatabricks → **user_impersonation** (admin consent recommended)           |
| Google Workspace    | BigQuery   | `https://www.googleapis.com/auth/bigquery`                                     |
| Other OIDC (Okta, …) | Databricks | Requires a [Databricks federation policy](https://docs.databricks.com/aws/en/admin/users-groups/best-practices-federation.html) on the workspace |

## Scopes

The connection's **Scopes** field is pre-filled from a per-connector,
per-flavour template. Defaults:

| Connector  | Flavour | Default scopes                                                           |
|------------|---------|--------------------------------------------------------------------------|
| BigQuery   | OIDC    | `https://www.googleapis.com/auth/bigquery`                               |
| BigQuery   | Azure   | `https://bigquery.googleapis.com/.default offline_access`                |
| Databricks | Azure   | `2ff814a6-3304-4ab8-85cb-cd0e6f879c1d/user_impersonation offline_access` |
| Databricks | OIDC    | `offline_access`                                                         |

Editing **Scopes** on a saved connection invalidates every user's stored
tokens — refresh tokens are bound to the original scope set. Each user sees
the consent prompt once more, then silent refresh resumes.

## Creating a connection

1. **Data → Databases**, right-click **BigQuery** or **Databricks**,
   **Add connection...**.
2. Fill in data-source-specific parameters (project, workspace URL, warehouse
   path, …).
3. Set **Auth** to **OAuth**.
4. Click **Test** — the consent popup opens. Sign in and consent.
5. Click **OK** to save.

End users need no setup. Each of them completes the same popup once on
first use.

## Security

- Tokens are stored encrypted in
  [credentials storage](../../../govern/access-control/access-control.md#credentials-management-system),
  scoped to each user's personal group — no cross-user leakage.
- The flow uses PKCE (RFC 7636) and a one-shot server-side `state` nonce
  (RFC 9700); returned scopes are compared against the requested set before
  storage.

## Troubleshooting

| Symptom                                        | Fix                                                                                                 |
|------------------------------------------------|-----------------------------------------------------------------------------------------------------|
| Popup: `redirect_uri_mismatch`                 | Add `https://<host>/oauth/connector/callback` to your IdP app's redirect URIs                       |
| Popup: `invalid_scope`                         | Grant the delegated permission in the IdP, or trim **Scopes** to match what's granted               |
| `consent_required` on every query              | Add `offline_access` to **Scopes** so a refresh token is issued                                     |
| Balloon returns after successful authorization | The IdP returned different scopes than requested — align the **Scopes** field with the IdP's policy |
| Databricks (OIDC): `401 Unauthorized`          | Configure a Databricks federation policy for your identity provider                                 |
| BigQuery: `Permission denied`                  | Grant the user `roles/bigquery.dataViewer` (or equivalent) in the Google Cloud project              |

## See also

- [BigQuery connector](bigquery.md)
- [Databricks connector](databricks.md)
- [Snowflake OAuth (JWT)](snowflake-oauth-jwt.md)
- [Data connection](../../access.md#data-connection)
