Testing of log in system comes down to checking system response to various combinations of entered data in **"Login or
E-mail"** and **"Password"** fields.

1. Open datagrok page in incognito mode.
2. Try the following combinations of data in fields when trying to log in.
   Should receive warning about failed log in and log in should not be performed:

* **"Login or E-mail"** - empty | **"Password"** - empty
* **"Login or E-mail"** - any data | **"Password"** - empty
* **"Login or E-mail"** - empty | **"Password"** - any data
* **"Login or E-mail"** - login that does not exist in system | **"Password"** - any data
* **"Login or E-mail"** - registered user login or email | **"Password"** - not correct for this login

3. Enter existing and correct combination of login and password.
   Successful login is performed. User enters.
   The *"Welcome"* view and has access to entities according to their permissions.

4. Log out by going to profile and click on ```Logout``` button.

5. Check the possibility to log in using account from third-party services - use **Login with Google**.

---
{
  "order": 1
}