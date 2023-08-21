---
title: "Configure authentication"
sidebar_position: 1
---

Datagrok supports [authentication](../../govern/authentication.md) using:

* [Login-password](#login-password-authentication)
* [Oauth](#oauth-authentication): 
  * [Google](#google-authentication)
  * [Facebook](#facebook-authentication)
  * [GitHub](#github-authentication)
* [OpenID](#openid-authentication)
* [LDAP](#ldap-authentication)

If you need coverage for a specific case, contact us on [community](https://community.datagrok.ai/) and we will discuss options.

## Login-password authentication

To configure login-password authentication for Datagrok:

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears.
2. To use the login-password method enable *Internal authentication* .
3. To disable signup, uncheck the *Signup Allowed* option. To disable OAuth signup, uncheck the *Allow Oauth Signup* option.
4. To restrict domains people can sign up from, use the *Signup Domains Whitelist* option. Separate several domains with commas.
5. To send welcome emails to new signed-up users enable the *Send Welcome Email* option.
6. To force users to use working emails enable the *Require Email Confirm* option.

To use login-password authentication, it is important to configure an [email service](configure-smtp.md) that delivers signup and password emails.

## OAuth authentication

Datagrok supports [Google](#google-authentication), [Facebook](#facebook-authentication) and [GitHub](#github-authentication) OAuth authentication.

### Google authentication

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears
2. To use the Google Oauth method enable *Google authentication* 
3. Set *Google Client Id* obtained from [Google](https://developers.google.com/identity/oauth2/web/guides/get-google-api-clientid)
4. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

### Facebook authentication

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears
2. To use the Facebook Oauth method enable *Facebook authentication* 
3. Set *Facebook Client Id* and *Facebook Secret*  obtained from [Facebook](https://help.vtex.com/tutorial/adding-a-client-id-and-a-client-secret-to-log-in-with-facebook--3R7rzXWG1GswWOIkYyy8SO)
4. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

### Github authentication

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears
2. To use the GitHub Oauth method enable *GitHub authentication* 
3. Set *GitHub Client Id* and *GitHub Secret* obtained from [GiHub](https://episyche.com/blog/how-to-create-oauth-client-id-and-client-secret-for-github)
4. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

## OpenID authentication

Datagrok supports the OpenID protocol to allow users to be authenticated using OpenID providers, for example, Azure AD.

To setup OpenID authentication:

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears.
2. To use the OpenID method enable *Open Id authentication* 
3. Get a well-known-configuration route and set it to **Open Id Config Endpoint**. It should look
   like `https://login.datagrok.ai/.well-known/openid-configuration`
4. Set *Open Id Client Id* and *Open Id Secret* as in your OpenId provider
5. Set the *Open Id Code Challenge method* if you enabled authorization code encryption. In most cases, it is `S256`
6. To provide optional claims for the application set *Open Id Login Claim*, *Open Id Email Claim*, *Open Id First Name Claim*, and *Open Id Last Name Claim* 
7. To enable enable/disable OpenID auto-login use the *Open Id Auto Login* option
8. Make sure the correct **Web Root** is set in the left **Sidebar** > **Settings** > **Admin** section

## LDAP authentication

To enable LDAP or Active Directory authentication:

1. On the left **Sidebar** go to **Settings** > **Users and Sessions**. Authentication settings page appears.
2. Enable *Domain authentication* to use the LDAP method
3. Enable *Domain signup* to enable all users present on a domain controller to authenticate in the Datagrok platform.
   If the option is disabled, it is required to create the user in the Datagrok platform first to allow the user to log into the platform
4. Configure LDAP server address/DNS name in *LDAP Host*
5. Set *LDAP port*
6. Enable *LDAP SSL* if you use LDAPS on your server
7. Set *LDAP Base DN*. The Base DN is a starting point within the directory's [hierarchical structure](https://docs.oracle.com/cd/E19182-01/820-6573/ghusi/index.html) from which LDAP searches and operations are performed. For example `dc=datagrok,dc=ai`.
8. To set user and groups set *LDAP User DN*. For example `CN=USER-DATAGROK,OU=users,DC=datagrok,DC=ai`
9. Set *LDAP User password*
