---
title: "Configure authentication"
sidebar_position: 1
---

Datagrok supports many [authentication](../../govern/access-control/access-control.md#authentication) methods, including popular methods such as SSO and OAuth:

* [General](#general-login-password-authentication) (login-password)
* [LDAP](#ldap-authentication)
* [Google](#oauth-authentication)
* [OpenID](#openid-authentication)
* [SAML](#saml-authentication)
* [IAP](#iap-authentication)
* GitHub
* Facebook

You can enable all methods separately or combined.

<!-- markdownlint-disable no-bare-urls -->
If supported authentication methods do not work for you, contact us on info@datagrok.ai, and we will discuss options for your specific case.
<!-- markdownlint-enable no-bare-urls -->

## General (login-password) authentication

General (login-password) authentication is the most basic method to authenticate users with Datagrok.

To configure login-password authentication:

1. Go to the **Settings** > **Users and Sessions**. This section contains all authentication settings.
2. To use the login-password method, enable '_Internal authentication_' in General section
3. To disable signup uncheck '_Signup Allowed_' option
4. To restrict from which domains people can sign up to the platform, use the '_Signup Domains Whitelist_' option. You can set several domains separated with commas.
5. To force people to use active emails, enable the '_Require Email Confirm_' option .

For login-password authentication, it is important to [configure an email service](configure-smtp.md) that will deliver signup, welcome, confirmation and forgot password emails.

### Add users

To create [user](../../govern/access-control/users-and-groups#users):

1. On the **Sidebar** go to **Manage** > **Users**. 
2. On the **Toolbox** click **Add User**. Create new user dialog appears. 
3. Fill all input fields and click **OK**. New User profile appears. Click **Save** on the **Top Bar**.

Use [user groups](../../govern/access-control/users-and-groups#groups) to manage user permissions inside platform.

## LDAP authentication

Datagrok integrates with your LDAP or Active Directory server enabling the smooth domain authentication mechanism across all your services.

1. Go to the **Settings** > **Users and Sessions**. This section contains all authentication settings.
2. To use the LDAP method, enable '_Domain authentication_' 
3. Enable 'Domain signup' to enable all users present on a domain controller to authenticate in the Datagrok platform.
   If the option is disabled, it is required to create the user in the Datagrok platform first to allow the user to log
   into the platform
4. Configure LDAP server address/DNS name
5. Set LDAP server port
6. Enable LDAP SSL if you use LDAPS on your server
7. Set LDAP Base DN. It should look like `dc=datagrok,dc=ai`.
8. Set LDAP User DN. It should look like `CN=USER-DATAGROK,OU=users,DC=datagrok,DC=ai`
9. Set LDAP User password

>Note: To ensure only domain-managed users can access the platform:
>
>1. Disable 'Signup Allowed' to prevent unauthorized users from registering directly on Datagrok.
>2. Enable 'Signup Enabled' in _'Domain Authentication'_ to allow new users already registered in the organization's LDAP or Active Directory (AD) system to log in.

## Oauth authentication

Datagrok supports Google, Facebook, and GitHub OAUTH authentication.

1. Go to the Datagrok Settings section 'Users and Sessions'; this section contains all authentication settings.
2. Enable 'Google authentication' to use the Google Oauth method (or another provider)
3. Set 'Client Id' and 'Secret' if applicable. You can get it from your OpenID provider
4. Make sure the correct Web Root is set in 'Admin' section

Consider using a general OpenID authentication as it is more flexible.

## OpenID authentication

This is the preferred and most powerful way to integrate with an external identity provider.  
Datagrok supports the OpenID protocol to allow users to be authenticated using OpenID providers, for example, Azure AD.

1. Go to the Datagrok Settings section `Users and Sessions`; this section contains all authentication settings.
2. Enable `Open Id authentication` to use the OpenID method
3. Get a well-known-configuration route and set it to 'Open Id Config Endpoint'. It should look
   like `https://login.datagrok.ai/.well-known/openid-configuration`
4. Set `Client Id` and `Client Secret` as in your OpenId provider. 
Datagrok supports a certificate-based authentication method for Azure AD. If you prefer it against a plain secret, 
generate a certificate, sign it with a private RSA key, and upload it to the Datagrok. Then upload the certificate to your Azure AD application settings. 
5. Set the `Code Challenge method` if you enabled authorization code encryption. In most cases, it is `S256`
6. Enable `Auto Login` option to forward users to authentication automatically without showing the login form.
7. Make sure the correct Web Root is set in `Admin` section
8. Enable `Keep Token` mode if you want to enable seamless integration with other services. 
Datagrok will request `offline_access` scope and keep encrypted external token in session metadata. 

### Keep Token

Default authentication token expiration is 1 hour for OpenID, and it only can be used for validating user identity.

To be able to use the token for accessing external services, enable `Keep Token` in OpenID settings.

When `Keep Token` is enabled, Datagrok requests offline access from the OpenID provider during authentication. 
This allows Datagrok to securely get a refresh token in addition to the access token.

Datagrok automatically refreshes the external OpenID token together with its own internal session token.
Datagrok keeps a refresh token in the browser storage.
Token refresh happens every 10 minutes and when a user session starts, ensuring uninterrupted user sessions and seamless integration with external services.

Datagrok applies different handling strategies depending on the token type, ensuring that server-side actions always require explicit user intent.
If the external token is a JWT, the token signature is embedded into the Datagrok user token on the client side.
The JWT payload is encrypted and stored in the database.
Datagrok cannot perform external actions autonomously without token signature.

If the external token is opaque, the token is stored only inside the user JWT on the client side.
The server never has direct access to the token contents.

That means, external tokens are never usable by the server without a user-initiated request

Long-lived sessions are supported without exposing privileged credentials, and Datagrok can integrate seamlessly with
external systems while preserving user-controlled authorization boundaries.

Right now user OpenID authentication is supported by [BigQuery](../../access/databases/connectors/bigquery.md) and [Databricks](../../access/databases/connectors/databricks.md) providers.

## SAML authentication

1. Go to Datagrok Settings section 'Users and Sessions'
2. Enable SAML authentication
3. Copy ACS URI and Entity ID to SAML provider (i.e. Google or Azure AD)
4. Copy SSO URI to Datagrok as IdP Endpoint and Certificate. Make sure you switch to a multiline edit mode when copying certificate.
5. Make sure you have mappings in your authentication provided: email, first_name, last_name 

## IAP authentication

Datagrok supports [Google Identity-Aware Proxy (IAP)](https://cloud.google.com/security/products/iap) out of the box. 
Configure Identity-Aware Proxy for Datagrok server for automatic login.

Datagrok automatically detects `x-goog-iap-jwt-assertion` header, validates the token using Google keys, and authenticates user.