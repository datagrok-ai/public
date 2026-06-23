---
feature: general
target_layer: manual-only
coverage_type: regression
produced_from: split
original_path: public/packages/UsageAnalysis/files/TestTrack/General/profile-settings.md
split_date: 2026-06-16
related_bugs: []
manual_only_reason: |
  Profile-photo upload goes through the Dart-only path user.setPicture() (image
  decode + resize in Dart) with no public JS API, and the "Change password..."
  modal's mismatch validation is a client-side Dart Modal with no name= hooks —
  both are fragile to drive and not exposed to the JS API used by Playwright
  specs. The name-edit (which calls dapi.users.save) is automated in
  profile-settings-spec.ts; the photo and password paths remain manual.
---

# Profile Settings — photo & password (manual)

Manual companion to `profile-settings.md`. The profile **name edit** and its
persistence are automated by `profile-settings-spec.ts`. This file covers the
photo-upload and password-validation paths, which are Dart/UI-only.

## Preconditions

- Logged in. Open the profile view (avatar → your name, or `/u/<login>`).
- Have ready: a few valid image files in different formats (e.g. `.png`,
  `.jpg`, `.gif`) and one non-image file (e.g. a `.txt` or `.csv`) for the
  negative case.

## Steps

### A. Profile photo

1. Click the profile picture area to change it; upload a valid image (positive
   scenario). Repeat with different popular image formats — each is accepted and
   the avatar updates.
2. Attempt to upload a file with an **incorrect / non-image format**. The photo
   should **not** change (no broken avatar; the invalid file is rejected).

### B. Change password (negative validation only)

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

---
{
  "order": 7
}
