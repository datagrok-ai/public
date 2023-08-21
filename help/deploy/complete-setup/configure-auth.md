---
title: "Configure authentication"
sidebar_position: 1
---

Datagrok supports authentication using:

* [Login-password](#login-password-authentication)
* [Oauth](#oauth-authentication): 
  * [Google](#google-authentication)
  * [Facebook](#facebook-authentication)
  * [GitHub](#github-authentication)
* [OpenID](#openid-authentication)
* [LDAP](#ldap-authentication)

You can enable all authentication methods separately or combined.

If supported authentication methods do not work for you, contact us on <info@datagrok.ai> or [our community](https://community.datagrok.ai/) and we will discuss options for your specific case.

[Learn more about authentification](../../govern/authentication.md).

## Common options

You can display a text on the login window, like greetings words etc.
To do so set the *Prompt* option in **Settings** > **Users and Sessions** > **General**:
<p align="center">
  <img alt="login-prompt" src={require('./login-prompt.png').default} width="300px"/>
</p>

## Login-password authentication

To configure login-password authentication for Datagrok:

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears.
2. Enable *Internal authentication* to use the login-password method.
3. Adjust the *Signup* option. You can disable or enable it using the *Signup Allowed* and *Allow Oauth Signup* checkboxes.
4. For better security you can restrict from which domains people can sign up to the platform. To do so use the *Signup Domains Whitelist* option. You can
   set several domains separated with commas
5. Enable the *Send Welcome Email* option to send welcome emails to new signed-up users.
6. Enable the *Require Email Confirm* option to force users to use working emails.

To use login-password authentication, it is important to configure an [email service](configure-smtp.md) that delivers signup and password emails.

## OAuth authentication

Datagrok supports [Google](#google-authentication), [Facebook](#facebook-authentication) and [GitHub](#github-authentication) OAuth authentication.

### Google authentication

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears
2. Enable *Google authentication* to use the Google Oauth method
3. Set *Google Client Id*. You can get it from [Google](https://developers.google.com/identity/oauth2/web/guides/get-google-api-clientid)
4. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

### Facebook authentication

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears
2. Enable *Facebook authentication* to use the Facebook Oauth method
3. Set *Facebook Client Id* and *Facebook Secret* if applicable. You can get it from [Facebook](https://help.vtex.com/tutorial/adding-a-client-id-and-a-client-secret-to-log-in-with-facebook--3R7rzXWG1GswWOIkYyy8SO)
4. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

### Github authentication

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears
2. Enable *GitHub authentication* to use the Facebook Oauth method
3. Set *GitHub Client Id* and *GitHub Secret* if applicable. You can get it from [GiHub](https://episyche.com/blog/how-to-create-oauth-client-id-and-client-secret-for-github)
4. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

## OpenID authentication

Datagrok supports the OpenID protocol to allow users to be authenticated using OpenID providers, for example, Azure AD.

To setup OpenID authentication:

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears.
2. Enable *Open Id authentication* to use the OpenID method
3. Get a well-known-configuration route and set it to **Open Id Config Endpoint**. It should look
   like `https://login.datagrok.ai/.well-known/openid-configuration`
4. Set *Open Id Client Id* and *Open Id Secret* as in your OpenId provider
5. Set the *Open Id Code Challenge method* if you enabled authorization code encryption. In most cases, it is `S256`
6. Set *Open Id Login Claim*, *Open Id Email Claim*, *Open Id First Name Claim*, and *Open Id Last Name Claim* to provide optional claims for the application
7. You can enable/disable OpenID auto-login using the *Open Id Auto Login* checkbox
8. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

## LDAP authentication

You can integrate Datagrok with your LDAP server by enabling the smooth domain authentication mechanism across all your
services. You can also use Active Directory to authenticate by this method.

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears.
2. Enable *Domain authentication* to use the LDAP method
3. Enable *Domain signup* to enable all users present on a domain controller to authenticate in the Datagrok platform.
   If the option is disabled, it is required to create the user in the Datagrok platform first to allow the user to log into the platform
4. Configure LDAP server address/DNS name
5. Set LDAP server port
6. Enable LDAP SSL if you use LDAPS on your server
7. Set LDAP Base DN. It should look like `dc=datagrok,dc=ai`.
8. Set LDAP User DN. It should look like `CN=USER-DATAGROK,OU=users,DC=datagrok,DC=ai`
9. Set LDAP User password
