---
feature: general
realizes_atlas: []
realizes: []
priority: p2
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  Profile-photo upload goes through the Dart-only path user.setPicture() (image
  decode + resize in Dart) with no public JS API, and the "Change password..."
  modal's mismatch validation is a client-side Dart Modal with no name= hooks —
  both are fragile to drive and not exposed to the JS API used by Playwright
  specs. The name-edit (which calls dapi.users.save) is automated in
  profile-settings-spec.ts; the photo and password paths remain manual.
related_bugs: []
---

# Profile Settings — photo & password (manual)

Manual companion to `profile-settings.md`. The profile **name edit** and its
persistence are covered by the automated companion. This file covers the
photo-upload and password-validation paths.

## Preconditions

- Logged in. Open the profile view (avatar → your name, or `/u/<login>`).
- Have ready: a few valid image files in different formats (e.g. `.png`,
  `.jpg`, `.gif`) and one non-image file (e.g. a `.txt` or `.csv`) for the
  negative case.

## Steps

### Profile photo

1. Click the profile picture area to change it; upload a valid image (positive
   scenario). Repeat with different popular image formats — each is accepted and
   the avatar updates.
2. Attempt to upload a file with an **incorrect / non-image format**. The photo
   should **not** change (no broken avatar; the invalid file is rejected).

### Change password (negative validation only)

1. Open **Change password...**.
2. Enter a **new password** and a **non-matching** confirmation → on OK, a
   "Password confirmation doesn't match" message appears and the password is
   **not** changed.
3. Enter a **wrong current password** with matching new/confirm values → on OK,
   a "Failed to change password" message appears and the password is **not**
   changed.

> Note: do not perform a *successful* password change against a shared test
> account — it would lock out other automated tests that authenticate as that
> user.

## Cleanup

- Restore the original avatar: re-upload the profile photo the account had
  before step A1 (the negative A2 upload leaves it unchanged).
- The password needs no revert — both B scenarios are negative cases and must
  leave it unchanged; if a change slipped through by mistake, change it back
  to the original immediately.

---
{
  "order": 7,
  "datasets": []
}
