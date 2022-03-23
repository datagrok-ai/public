<!-- TITLE: Tests: Login -->
<!-- SUBTITLE: -->

# Tests: Login

Platform provides multi-user system. Users can create their own accounts or use third-party services accounts (Google,
GitHub, Facebook) for login

## Testing scenario

Testing of log in system comes down to checking system response to various combinations of entered data in **"Login or
E-mail"** and **"Password"** fields. With the following combinations of data in fields when trying to log in, user
should receive warning about failed log in and log in should not be performed:

* **"Login or E-mail"** - empty | **"Password"** - empty
* **"Login or E-mail"** - any data | **"Password"** - empty
* **"Login or E-mail"** - empty | **"Password"** - any data
* **"Login or E-mail"** - login that does not exist in system | **"Password"** - any data
* **"Login or E-mail"** - registered user login or email | **"Password"** - not correct for this login

When user enter existing and correct combination of login and password, successful login is performed. User enters
the *"Welcome"* view and has access to entities according to their permissions.

To log out, user must go to his profile and click on ```Logout``` button.

User can use existing email in combination with correct password to log in.

Also need testing possibility of logging using accounts from third-party services. In this case, if user agrees to use
his Google (Facebook, GitHub) account to log in to platform, then he can be logged in.

See also:

* [User](../../govern/user.md)
* [Authentication](../../govern/authentication.md)
* [Authorization](../../govern/authorization.md)
* [Login Auto Test](../selenium/login-test.side)
