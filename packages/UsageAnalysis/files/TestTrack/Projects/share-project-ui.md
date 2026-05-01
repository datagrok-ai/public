---
feature: projects
companion_to: share-project.md
companion_spec: share-project-spec.ts
ui_only: true
moved_steps_from_canonical: ["Step 4 sub-bullet (email)", "Step 5", "Step 7", "Step 8"]
generated_at: 2026-05-01
generated_by: D3 bucket-b classification (see decision-log + coverage-gaps/projects.md)
---

# Share Project — UI-only manual companion

This file captures the UI-only verification steps from
`share-project.md` that cannot be cleanly automated via JS API. Step
numbering preserved from the canonical `share-project.md` for
cross-reference.

These manual steps are NOT exercised by `share-project-spec.ts`. Human
QA must run them as a companion when validating the share scenario
end-to-end.

## Pre-requisite

Either: (a) `share-project-spec.ts` has run successfully to create the
fixture project (name pattern `AutoTest-Share-<timestamp>`), OR (b)
the saved `demog` project from `upload-project.md` exists on the test
environment.

## Manual verification steps

### Step 4 (email sub-bullet) — Share via email to an unregistered email recipient

1. With the saved test project selected (Browse > Dashboards), if not
   already in the Share dialog: right-click the project → **Share**.
2. In the Share dialog, enter an email address that does NOT correspond
   to any existing Datagrok account at run time (e.g.
   `qa-invite-<timestamp>@example.com`).
3. Choose access level (View-and-Use / Edit / Full).
4. Confirm/Submit the share.

**Why UI-only:** the auto-create-user-account-from-email behavior is
a side-effect of the Share dialog's email-handling logic. JS API path
(`grok.dapi.permissions.grant(p, '<email>', false)`) does not reliably
trigger the account-provisioning logic on current dev (returns
`NoSuchMethodError: getter 'id' was called on null` — see Wave 1a
share-project-spec.verdict.yaml `validator_rounds`).

### Step 5 — Verify on Context Panel — Sharing tab (auto-created user listed)

1. After the share submits, navigate to the project's Context Panel >
   **Sharing** tab.
2. Verify a NEW user account has been auto-created for the
   unregistered email address — visible in the Sharing tab recipient
   list (display name may be the email itself or a derived handle).
3. Verify the recipient's permission level matches what was selected
   in Step 4.

**Why UI-only:** confirming that a NEW user entity was
auto-provisioned (vs an existing user already in the system) is most
reliably done by visually inspecting the Sharing tab AND cross-checking
that the recipient's login was absent from the pre-share user list.
While `dapi.permissions.list` returns share relations, it does not
distinguish "new user just created" from "existing user".

### Step 7 — Right-click the project and select Details

1. Browse > Dashboards. Locate the test project.
2. Right-click the project node.
3. From the context menu, select **Details**.

**Why UI-only:** context menu items have no `name=` attributes
(`grok-browser/references/projects.md:242`). Locating "Details" by
visible text is fragile; the right-click → Details dispatch flow is
best exercised manually when verifying Context Panel render quality
(Step 8).

### Step 8 — Verify on Context Panel — all tabs render correctly

After Step 7, the project's Details / Context Panel is visible. Review
**every tab** on the Context Panel:

- **Sharing** — recipients list renders with avatars/icons; permission
  level labels visible; empty-state placeholder if no recipients.
- **History** — audit entries chronological; timestamps formatted;
  entries readable.
- **Properties** — project metadata key-value pairs aligned; no
  blank/error values.
- **Activity** / **Chats** / **Tables** (if present) — each renders
  appropriately for content presence/absence.

For each tab, check (subjective render quality):
- Layout: no overlapping elements, no truncated text, no missing labels.
- Typography: consistent font sizing/weight; alignment matches
  surrounding panels.
- No raw error messages, broken icons, or rendering glitches.

**Why UI-only:** "render correctly" is a subjective render-quality
assessment — layout, typography, no-glitch criterion is human judgment
without an explicit per-tab automated criterion.

## Cleanup responsibility

After the manual run, clean up the auto-created user account from
Step 4 to keep the environment idempotent:

- Via `grok s users delete <login>` CLI, OR
- Via Manage > Users in the Datagrok UI as an admin.

Log the cleanup in `share-project-run.md` under the dated entry for
this manual run.

## Sign-off

PASS — all 4 manual steps completed; no visual regressions vs prior
dev release.
FAIL — list affected step(s) + screenshot; log in `share-project-run.md`
under a new dated entry.
