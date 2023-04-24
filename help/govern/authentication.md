---
title: "Authentication"
---

Out-of-the-box, Datagrok offers the following authentication methods:

* [Login-password](#login-password)
* OAUTH (Google, Facebook, GitHub)
* Single Sign-on
* Active Directory

Enterprise customers might prefer to use a custom SSO (single sign-on) scheme. We can accommodate these needs by
developing a customer-specific integration.

## Login-password

Datagrok uses login and password to authenticate users. We keep your password salted with random data and encrypted with
1024xSHA-256 algorithm, so it cannot be read from our system. If you forgot your password - the only way to get access
to Datagrok is to reset your password, using link on login form.

## Authentication details

When a user enters logs in, the login and password pair is passed to the server, and, if the password hash matches the stored
hash, a session token is generated. Now, every API call should be made with http header "Authorization: token",
where "token" is the session token. After logging out this token won't work anymore.

Datagrok does not keep your password anywhere, after your log in.

![Authentication UML Diagram](../uploads/features/login-signup.png "Authentication UML Diagram")
[Authentication UML Diagram draw.io](https://github.com/datagrok-ai/public/tree/master/help/uploads/features/login-signup.drawio)

See also:

* [Authorization](authorization.md)
* [Security](security.md)
* [Configure authentication](../develop/admin/deploy/configure-auth.md)
