<!-- TITLE: Tests: Registration -->
<!-- SUBTITLE: -->

# Tests: Registration

New user needs to fill information in the following fields of registration form: "Login", "First Name", "Last Name", "
New Password", "Retype Password". Admin can add new users.

## User signing up

1. Click on "Sign up" tab on start screen.

1. Click on "Sign Up" button without filling fields

* Warning "Failed to Sing up. Fill in required fields"

1. Enter incorrect Email address in appropriate field. (e.g. "test.mail.com", "test@mail", "
   test@.mail")

* Warning message about incorrect Email

1. Enter valid Email and then enter incorrect data in "Login" field (less than six symbols or with special symbols)

* Warning message about incorrect Login

1. Fill in "Password" and "Retype password" fields with different data

* Warning message about mismatched passwords

1. Use already existing email

* "@email already exists" warning

1. Fill in all fields with valid data and click "Sign Up"

* New user added (to check: ```Admin | Users```)

1. Check email address provided during registration

* "Welcome" mail received

## Administrator creating a new user account

Pre-requisites: an admin has to have the privilege to add new users.

1. Open ```Users``` from ```Admin``` menu.

1. Click on ```Add user``` button

* "Create new user" dialog is open

1. Click ```OK``` without filling fields

* Warning "Failed to Sing up. Fill "E-mail" field"

1. Enter incorrect Email address in appropriate field. (e.g. "test.mail.com", "test@mail", "
   test@.mail")

* Warning message about incorrect Email

1. Enter valid Email and then click ```OK```

* New user setup page is open
* Since, "Login" field was not filled, login corresponds to first part of email

1. Click on "groups". Select group to add new user to it. (there must be corresponding global privileges)

1. Click on ```SAVE``` button

* New user added
* New user enters selected groups
* By using specified data, user can log in

See also:

* [User](../../govern/user.md)
* [Authentication](../../govern/authentication.md)
* [Authorization](../../govern/authorization.md)
