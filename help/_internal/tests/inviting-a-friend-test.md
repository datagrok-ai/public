<!-- TITLE: Tests: Inviting a friend -->
<!-- SUBTITLE: -->

# Tests: Inviting a friend

The system allows you to invite a friend by email.

## Testing scenarios

1. Enter an incorrect Email address in the Invite dialog. (e.g. "test.mail.com", "test@mail", "
   test@.mail")

* Invitation not sent. Warning message about incorrect Email.

1. Send an invitation to a valid email address

* Invitation sent to email successful. Invited user received a letter with invitation link. After following a link, a
  new user was created.

See also:

* [Registration test](../tests/registration-test.md)
