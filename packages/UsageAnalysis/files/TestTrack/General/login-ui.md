---
feature: general
realizes_atlas: []
realizes: []
priority: p2
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  Third-party OAuth (Login with Google) redirects to an external Google-owned
  consent screen whose DOM and bot-detection are outside Datagrok's control and
  not stable to automate. The username/password positive+negative matrix and
  logout are already fully automated in General/login.test.ts; only the Google
  sign-in path remains manual.
related_bugs: []
---

# Login with Google (third-party auth)

Manual companion to `login.md`. The local-credential login matrix (empty /
wrong / boundary / special inputs, valid login, logout) is covered by the
automated companion. This file covers only the external OAuth path.

## Preconditions

- A browser session that is **logged out** (use a fresh incognito window).
- A Google account that is registered with / permitted to sign in to the target
  Datagrok instance.

## Steps

1. Open the Datagrok page in incognito mode — the login screen is shown.
2. Click **Login with Google**.
3. On the Google consent screen, choose the account and approve.
4. Expected result: the browser returns to Datagrok and the user is signed in —
   the *Welcome* view is shown and entities are accessible per the user's
   permissions.
5. Log out (profile → **Logout**) and confirm the login screen returns.

---
{
  "order": 1,
  "datasets": []
}
