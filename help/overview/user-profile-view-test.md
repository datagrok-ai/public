<!-- TITLE: Tests: Profile Settings -->
<!-- SUBTITLE: -->

# Tests: Profile Settings

[Users](../govern/user.md) can change the information in their profiles. Such as profile photo, password, first and last names.

## Testing scenarios

1. Change user photo using an incorrect file format
   * Display warning message about the wrong file format

1. Change user photo using an different formats of image files
   * All popular image formats are supported 

1. Edit name

1. Change password using negative scenario (invalid old password, mismatched passwords)
   * Password is not changed, the appropriate validation message should be shown on the dialog.

1. Use positive scenario to edit profile information (user photo, password, name)
   * Profile information is changed by entry data

See also:
  * [User](../govern/user.md)
  * [Registration test](../tests/registration-test.md)
