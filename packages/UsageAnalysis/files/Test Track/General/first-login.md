### First Login – Initial Workspace Validation

#### Preconditions:
- A new user account exists.User has never logged into Datagrok before.

#### Steps & Expected Results

1. Login
- Log in with the new account. Login succeeds without errors, redirects to the default start page.
2. Browse Tree – Tutorials Expansion
- Browse. Tutorials > Cheminformatics folder is expanded by default.
- Other tutorial categories remain collapsed.
3. Spaces Folder Check
- Navigate to Browse > Spaces. Only default/system spaces are visible. Currently Spaces folder should be empty. No other users’ content is present.
4. MyStuff Folder Check
- Navigate to MyStuff > My Files. Folder is empty, only one readme file is present.
5. Files > Dashboards. Each demo project opens successfully without errors, missing data.
6. UI / Landing Page Elements
- Check visible elements on landing page after first login. No errors, warnings, or broken widgets.

#### Postconditions:
- User remains logged in after tab reloading.
- No extraneous data or other users’ artifacts are visible.
- Demo/tutorial resources open without issues.

---
{
  "order": 10
}
