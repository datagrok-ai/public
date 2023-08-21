---
title: "Configure authentication"
sidebar_position: 1
---

Datagrok supports authentication using:

* [Login-password](#login-password-authentication)
* [Oauth](#oauth-authentication)
* [OpenID](#openid-authentication)
* [LDAP](#ldap-authentication)

You can enable all authentication methods separately or combined. To accommodate any enterprise needs, we can develop a customer-specific integration. 

[Learn more about authentification](../govern/authentication.md).

## Common options

To display something on the login window, set the **Prompt** option
<p align="center">
  <img alt="login-prompt" src={require('./login-prompt.png').default} width="300px"/>
</p>

Enable the **Send Welcome Email** option to send welcome emails to new signed-up users.

## Login-password authentication

To configure login-password authentication for Datagrok:

1. Go to the Datagrok Settings section **Users and Sessions**; this section contains all authentication settings.
2. Enable **Internal authentication** to use the login-password method
3. Signup option is alterable; you can disable it using the **Signup Allowed** option
4. To restrict from which domains people can sign up to the platform, use the **Signup Domains Whitelist** option. You can
   set several domains separated with commas
5. Enable the **Require Email Confirm** option to force people to use working emails.

To use login-password authentication, it is important to configure an [Email service](configure-smtp.md) that will deliver signup, welcome,
confirmation and forgot password emails.

## Oauth authentication

Datagrok supports **Google**, **Facebook** and **GitHub** OAUTH authentication.

1. Go to the Datagrok Settings section **Users and Sessions**; this section contains all authentication settings.
2. Enable **Google authentication** to use the Google Oauth method (or another provider)
3. Set **Client Id** and **Secret** if applicable. You can get it from your OpenID provider
4. Make sure the correct Web Root is set in **Admin** section

## OpenID authentication

Datagrok supports the OpenID protocol to allow users to be authenticated using OpenID providers, for example, Azure AD.

To setup OpenID authentication:

1. Go to the Datagrok Settings section **Users and Sessions**; this section contains all authentication settings.
2. Enable **Open Id authentication** to use the OpenID method
3. Get a well-known-configuration route and set it to **Open Id Config Endpoint**. It should look
   like `https://login.datagrok.ai/.well-known/openid-configuration`
4. Set **Open Id Client Id** and **Open Id Secret** as in your OpenId provider
5. Set the **Open Id Code Challenge method** if you enabled authorization code encryption. In most cases, it is `S256`
6. Set **Open Id Login Claim**, **Open Id Email Claim**, **Open Id First Name Claim**, and **Open Id Last Name Claim** to provide optional claims for the application
7. You can enable OpenID auto-login using the **Open Id Auto Login** option
8. Make sure the correct Web Root is set in **Admin** section

## LDAP authentication

You can integrate Datagrok with your LDAP server by enabling the smooth domain authentication mechanism across all your
services. You can also use Active Directory to authenticate by this method.

1. Go to the Datagrok Settings section **Users and Sessions**; this section contains all authentication settings.
2. Enable **Domain authentication** to use the LDAP method
3. Enable **Domain signup** to enable all users present on a domain controller to authenticate in the Datagrok platform.
   If the option is disabled, it is required to create the user in the Datagrok platform first to allow the user to log into the platform
4. Configure LDAP server address/DNS name
5. Set LDAP server port
6. Enable LDAP SSL if you use LDAPS on your server
7. Set LDAP Base DN. It should look like `dc=datagrok,dc=ai`.
8. Set LDAP User DN. It should look like `CN=USER-DATAGROK,OU=users,DC=datagrok,DC=ai`
9. Set LDAP User password
