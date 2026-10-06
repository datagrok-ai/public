---
feature: general
realizes_atlas: []
realizes: [views.space, views.files]
priority: p2
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Requires a never-logged-in account (one-shot state that an automated run
  would consume) and visual judgment of the first-login workspace: default
  tree expansion, empty folders, landing-page widgets.
related_bugs: []
---

# First login — initial workspace validation

## Preconditions

- A user account that has **never logged into Datagrok** on this stand.
  Create one for the run: Browse > Platform > Users > Add user (or ask an
  admin), e.g. `test_first_login_<date>`.

## Steps

1. Log in with the new account — login succeeds without errors and lands on
   the default start page
2. Open **Browse** — under Tutorials, the **Cheminformatics** folder is
   expanded by default; the other tutorial categories are collapsed
3. Go to **Browse > Spaces** — only default/system spaces are visible; no
   other users' content is present
4. Go to **My stuff > My Files** — the folder contains only the readme file,
   nothing else
5. Go to **Browse > Dashboards** and open the first two demo dashboards
   listed — each opens without errors or missing data
6. Return to the start page — the landing widgets render without error
   balloons or broken (blank/red) widget cards; empty states are acceptable
   for a new user; no errors in the browser console
7. Reload the browser tab — the user is still logged in, the platform
   restores without re-asking credentials

## Cleanup

Delete the test account created for the run (Browse > Platform > Users), or
hand it back to the admin who provisioned it. Note: the account is spent —
a second run needs a fresh never-logged-in account.

---
{
  "order": 2,
  "datasets": []
}
