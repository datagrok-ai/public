---
feature: users
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: []
realized_as:
  - users.test.ts
related_bugs: []
---

# Users management — manual test cases (dev / existing-users set)

**Instance:** https://dev.datagrok.ai/
**Object under test:** Users View (Browse > Platform > Users) — listing, NEW dialogs, search,
view modes, context menu, user profile, and admin management actions (block, edit groups/roles).
**Type:** UI / manual (staged for a Playwright port)
**Sources:** wiki `govern/access-control/users-and-groups.md`, `datagrok/navigation/views/user-profile-view.md`;
live recon on dev (see `user_groups_tests_context.md`).

---

## Two test sets — read first

These tests come in **two sets** that differ only in **user creation**:

- **This set (reddata, runs on dev):** never creates users. Users **cannot be hard-deleted**, so
  creating them on dev would pile up junk accounts. All management cases here operate on **existing**
  `opavlenko<NNN>` test users, modified by the admin account `opavlenko+playwright`. Create-user cases
  (Users-05, Users-07) are **deferred** — see below.
- **CI/CD set (`public/playwright-public/user groups/`, fresh DB):** adds the full create-user /
  service-user E2E (Users-05, Users-07 active). Built when porting to CI/CD.

**Existing target users:** any disposable `opavlenko<NNN>` account from the Users list (e.g.
`opavlenko45`, `opavlenko9`, `opavlenko656`). **Do not** target `opavlenko+playwright` (the runner)
or `opavlenko+pwsharing` (used by sharing tests). Always revert any change in cleanup.

---

## Conventions

- **ID:** `Users-<NN>` — Test Track compatible.
- **Priority:** P1 (critical path) … P3 (edge / cosmetic).
- **Type:** `smoke` / `functional` / `regression` / `negative`.
- Each case: **Goal → Type → Preconditions → Steps (action → expected UI reaction) → Final check →
  Postconditions/cleanup → Selectors**.

**Global preconditions (unless a case says otherwise):**
- Logged in to https://dev.datagrok.ai/ as the **admin** account (`DATAGROK_LOGIN`); app loaded.
- Console-error capture is on (`pageerror`, `console.error`, error balloons) for "no errors" cases.
- Browse panel open; Context Panel open (F4).

**What counts as an error:** uncaught `pageerror`; `console.error`; a Datagrok error balloon; a
Dart/JS stack surfaced in the UI or Context Panel.

---

## 1. Users View — open and layout

### Users-01 — Open Users View from the Browse tree (P1)
**Goal:** Users View opens from Browse > Platform > Users.
**Type:** smoke
**Preconditions:** global.
**Steps:**
1. In the Browse tree, expand **Platform**. → The Platform group expands and shows subnodes (Admin, Plugins, Functions, **Users**, **Groups**, Roles, …).
2. Click **Users**. → The Users View opens in the center; the URL becomes `…/users`.
3. Inspect the view. → A list of user entities is shown; the header has the **NEW** button, a search field, view-mode toggles, and a `shown / total` count.

**Final check:** active view type is `users`; URL contains `/users`; the user list and the count are visible; no errors.
**Postconditions/cleanup:** none.
**Selectors:** `treeNodeByPath(['Platform','Users'])`, `GALLERY_COUNTS`, `NEW_USER_BUTTON`, `GALLERY_SEARCH`.

---

### Users-02 — Users View toolbar controls are present (P2)
**Goal:** all expected header controls render.
**Type:** smoke
**Preconditions:** Users View open.
**Steps:**
1. Inspect the header. → Present: **NEW** combo button; **Search users by name or by #tags** input; view toggles **brief / card / grid**; **Sort list**; **Toggle filters**; **Refresh**; count `N / N`.

**Final check:** every listed control is visible and enabled.
**Postconditions/cleanup:** none.
**Selectors:** `NEW_USER_BUTTON`, `GALLERY_SEARCH`, `VIEW_TOGGLE_BRIEF/CARD/GRID`, `SORT_LIST`, `TOGGLE_FILTERS`, `[name="icon-sync"]`, `GALLERY_COUNTS`.

---

## 2. NEW — create dialogs (non-creating checks only on dev)

### Users-03 — NEW menu lists all create options (P1)
**Goal:** the NEW combo exposes the three create options.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Click **NEW**. → A dropdown opens with exactly: **User...**, **Service User...**, **Invite a Friend...**.

**Final check:** the three items are present in that order.
**Postconditions/cleanup:** close the menu (Esc) — nothing is created.
**Selectors:** `NEW_USER_BUTTON`, `.d4-menu-popup .d4-menu-item-label`.

---

### Users-04 — "Create new user" dialog fields and cancel (P1)
**Goal:** the New User dialog has the right fields and CANCEL is non-destructive.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Click **NEW > User...**. → The **Create new user** dialog opens.
2. Inspect fields. → Inputs: **Email**, **Login**, **First Name**, **Last Name**; buttons **OK**, **CANCEL**.
3. Click **CANCEL**. → The dialog closes; no user is created; the count is unchanged.

**Final check:** dialog title is `Create new user`; all four fields present; after CANCEL the list/count are unchanged.
**Postconditions/cleanup:** none. **Never** click OK in this set (would create an undeletable user).
**Selectors:** `NEW_USER_BUTTON`, `.d4-dialog`, `[name="input-host-Email"|"input-host-Login"|"input-host-First-Name"|"input-host-Last-Name"]`, `.ui-btn`.

---

### Users-05 — Create a new user end-to-end (CI/CD set only) (P1)
**Goal:** a user created via the dialog appears in the list.
**Type:** functional
**Status:** **CI/CD set only — NOT run on dev.** Users have no hard delete, so a real create pollutes
the instance. The full create-user E2E (enter Email/Login/First/Last → OK → user appears → verify by
search → block in cleanup) lives in `public/playwright-public/user groups/users_manual_tests.md` and
runs only on a fresh-DB CI/CD instance.
**On dev:** covered indirectly by Users-04 (dialog fields + CANCEL) and Users-06 (validation).

---

### Users-06 — New user validation: invalid/empty input (P2)
**Goal:** the dialog rejects an empty/invalid email.
**Type:** negative
**Preconditions:** Users View open.
**Steps:**
1. Click **NEW > User...**. → Dialog opens.
2. Leave **Email** empty (or enter `not-an-email`) and click **OK**. → The dialog does **not** create a user: OK is blocked or a validation message/error balloon appears; the dialog stays open or shows the error.

**Final check:** no user is created on invalid input; the user is informed (disabled OK or error message). The count is unchanged.
**Postconditions/cleanup:** CANCEL the dialog. (If, despite invalid input, a user were created, that is a bug to report — do not retry blindly.)
**Selectors:** dialog inputs, `.ui-btn`, `BALLOON_CONTAINER`.

---

### Users-07 — Create a service user (CI/CD set only) (P2)
**Goal:** the Service User dialog creates a service account.
**Type:** functional
**Status:** **CI/CD set only — NOT run on dev** (same no-hard-delete reason as Users-05). The dialog
shape (single **Login** field, OK/CANCEL) is verified non-destructively on dev within Users-08-style
open+cancel; the actual create runs only on CI/CD.

---

### Users-08 — "Invite a Friend" dialog (UI only) (P3)
**Goal:** the Invite dialog has an Email field; no email is actually sent.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Click **NEW > Invite a Friend...**. → The **Invite a Friend** dialog opens with an **Email** field and OK/CANCEL.
2. Inspect the field, then click **CANCEL**. → Dialog closes; nothing is sent.

**Final check:** dialog title `Invite a Friend`; Email field present; CANCEL is non-destructive.
**Postconditions/cleanup:** none. **Do not** click OK (sends a real email).
**Selectors:** `NEW_USER_BUTTON`, `.d4-dialog`, `[name="input-host-Email"]`.

---

## 3. Search, view modes, sort, filters

### Users-09 — Search users by name (P1)
**Goal:** typing a query filters the list and reflects in the URL.
**Type:** functional
**Preconditions:** Users View open; note total `N`.
**Steps:**
1. Type a known user fragment (e.g. `opavlenko`) in the search field. → The list narrows; the count drops below `N`; the URL gains `?q=opavlenko`.
2. Clear the search. → The full list returns; count back to `N / N`; `?q=` removed.

**Final check:** filtered count < total while searching; restored on clear; URL tracks the query.
**Postconditions/cleanup:** clear search.
**Selectors:** `GALLERY_SEARCH`, `GALLERY_COUNTS`, URL assertion `/users?q=…/`.

---

### Users-10 — Search users by #tag (P3)
**Goal:** `#tag` search syntax is accepted.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Type `#` in the search field. → A tag suggestion/autocomplete may appear.
2. Enter a tag query. → The list filters to users matching the tag (possibly empty); no errors.

**Final check:** tag search runs without error; the count updates.
**Postconditions/cleanup:** clear search.
**Selectors:** `GALLERY_SEARCH`, `GALLERY_COUNTS`.

---

### Users-11 — Switch list view modes (P2)
**Goal:** brief / card / grid toggles change the layout.
**Type:** functional
**Preconditions:** Users View open (brief by default).
**Steps:**
1. Click **Switch to card view**. → The list re-renders as cards; the card toggle becomes current.
2. Click **Switch to grid view**. → The list re-renders as a data grid (tabular).
3. Click **Switch to brief view**. → Back to the brief list.

**Final check:** each toggle visibly changes the layout and marks itself current (`d4-current`); no errors.
**Postconditions/cleanup:** return to brief view.
**Selectors:** `VIEW_TOGGLE_BRIEF/CARD/GRID` (`d4-current` class on active).

---

### Users-12 — Sort the user list (P3)
**Goal:** the Sort control reorders the list.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Click **Sort list**. → A sort options popup opens (e.g. by name / created).
2. Pick a sort option. → The list reorders accordingly; no errors.

**Final check:** the order changes after applying a sort; no errors.
**Postconditions/cleanup:** none.
**Selectors:** `SORT_LIST`, `.d4-menu-popup`.

---

### Users-13 — Toggle the filters panel (P2)
**Goal:** the filter toggle opens a filter panel for user properties.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Click **Toggle filters**. → A filter panel appears (e.g. by group / created / login).
2. Apply one filter. → The list narrows; the count drops; no errors.
3. Toggle filters off / reset. → The list restores.

**Final check:** filter panel opens; applying a filter changes the count; resetting restores it; no errors.
**Postconditions/cleanup:** reset and close filters.
**Selectors:** `TOGGLE_FILTERS`, filter panel, `GALLERY_COUNTS`.

---

## 4. User context menu and profile

### Users-14 — User context menu items (P1)
**Goal:** right-clicking a user shows the expected actions.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Right-click an existing user (e.g. `opavlenko45`). → A context menu opens with: **Details**, **Chat**, **Block**, **Groups...**, **Roles...**, **ID** (copy id), **Add to favorites**.

**Final check:** all listed items are present.
**Postconditions/cleanup:** close the menu (Esc) — no action executed.
**Selectors:** user card, `.d4-menu-popup .d4-menu-item`, `d4-name="Details|Chat|Block|Groups...|Roles...|Add to favorites"`.

---

### Users-15 — Open a user profile by double-click (P1)
**Goal:** double-clicking a user opens its profile view.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Double-click an existing user (e.g. `opavlenko45`). → The **user profile** view opens (view type `user_edit`), titled with the user's name.
2. Inspect the Context Panel. → Info panes are shown: **Personal**, **Roles**, **Member of**, **Projects**, **Activity**, **Privileges**, **Chats**, **Sticky meta**.

**Final check:** the profile view opens; the listed info panes are present; no errors.
**Postconditions/cleanup:** close the profile view.
**Selectors:** user card (dblclick), view type `user_edit`, `.grok-prop-panel .d4-accordion-pane[name="pane-Personal|pane-Roles|pane-Member-of|…"]`.

---

### Users-16 — Open a user profile via "Details" (P2)
**Goal:** the Details context action opens the same profile view.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Right-click an existing user → **Details**. → The user profile view opens (as in Users-15).

**Final check:** profile view opens; no errors.
**Postconditions/cleanup:** close the view.
**Selectors:** context menu `d4-name="Details"`, view type `user_edit`.

---

### Users-17 — User context-panel info panes (selecting in list) (P2)
**Goal:** clicking a user once populates the Context Panel with user info panes.
**Type:** smoke
**Preconditions:** Users View open; Context Panel open.
**Steps:**
1. Single-click an existing user. → The Context Panel header shows the user; panes appear: **Personal**, **Roles**, **Member of**, **Projects**, **Activity**, **Privileges**, **Chats**, **Sticky meta**.
2. Expand **Personal**. → Name, email, login, join date are shown.

**Final check:** the expected panes render; Personal shows identity fields; no errors.
**Postconditions/cleanup:** none.
**Selectors:** user card (click), `.grok-prop-panel .d4-accordion-pane[name="pane-…"]`.

---

## 5. Admin management actions — on existing users (reverted in cleanup)

### Users-18 — Edit an existing user's group memberships (Groups...) (P1)
**Goal:** the Groups... dialog adds/removes the user's group memberships.
**Type:** functional
**Preconditions:** Users View open; admin; a target existing user `opavlenko<NNN>`; a known target group exists (e.g. a `qa_autotest` group created by the Groups suite, or any non-personal group the user is not in). Record the user's current memberships first.
**Steps:**
1. Right-click the target user → **Groups...**. → The `<user> groups` dialog opens with a checkbox list of groups and SAVE/CANCEL.
2. Check a group the user is **not** currently in. → The checkbox becomes selected.
3. Click **SAVE**. → Dialog closes; the membership is saved.
4. Single-click the user and expand **Member of** in the Context Panel. → The newly added group is listed.

**Final check:** the user gains membership in the chosen group (visible in **Member of**).
**Postconditions/cleanup:** reopen **Groups...**, uncheck the group, SAVE — restore the original memberships.
**Selectors:** context menu `d4-name="Groups..."`, `.d4-dialog`, group checkboxes, `pane-Member-of`.

---

### Users-19 — Edit an existing user's roles (Roles...) (P2)
**Goal:** the Roles... dialog assigns/removes roles.
**Type:** functional
**Preconditions:** Users View open; admin; a target existing user `opavlenko<NNN>`. Record current roles first.
**Steps:**
1. Right-click the target user → **Roles...**. → The `<user> roles` dialog opens listing roles (each with a **Can assign** option) and SAVE/CANCEL.
2. Toggle a role on that the user does not have. → The role is selected.
3. Click **SAVE**. → Dialog closes; the role assignment is saved.
4. Single-click the user and expand **Roles** in the Context Panel. → The assigned role is shown.

**Final check:** the role appears in the user's **Roles** pane.
**Postconditions/cleanup:** reopen **Roles...**, toggle the role off, SAVE — restore the original roles.
**Selectors:** context menu `d4-name="Roles..."`, `.d4-dialog`, role toggles, `pane-Roles`.

---

### Users-20 — Block a user (and unblock it back) (P1)
**Goal:** Block disables an account; the account is then restored.
**Type:** functional
**Preconditions:** Users View open; admin; a **disposable** existing user `opavlenko<NNN>` (never `opavlenko+playwright` or `opavlenko+pwsharing`).
**Steps:**
1. Right-click the target user → **Block**. → An **"Are you sure? Block user …"** confirmation appears (buttons **YES** / **CANCEL**).
2. Click **YES**. → The user's server status becomes `blocked` (the block is async — verify by polling the status, e.g. `GET /api/users/{id}`.status).

**Final check:** the user's status is `blocked` after confirming.
**Observed on dev (important):** the Users gallery has **no UI "Unblock"** — the context menu does **not** flip to Unblock after blocking, and the user-profile view has no block/unblock control either. Unblocking is only possible via the server API (`POST /api/users/unblock` with the user JSON). So:
**Postconditions/cleanup:** **unblock via the API** (`POST /api/users/unblock`) in `finally`, and assert the status returns to `active`. (This is the documented exception to UI-first — there is no UI path to unblock.)
**Selectors:** context menu `d4-name="Block"`, confirm `button[name="button-YES"]`; status/unblock via API.

---

### Users-21 — Add an existing user to favorites (P3)
**Goal:** a user can be starred and shows under Favorites.
**Type:** functional
**Preconditions:** Users View open.
**Steps:**
1. Right-click an existing user → **Add to favorites**. → The user is starred (the Context Panel star becomes active when the user is selected).
2. Open **Browse > My stuff > Favorites** (or the Favorites sidebar). → The user appears in the favorites list.

**Final check:** the user is present in Favorites.
**Postconditions/cleanup:** remove the user from favorites (toggle the star off).
**Selectors:** context menu `d4-name="Add to favorites"`, `CONTEXT_PANEL_STAR`, Favorites list.

---

## 6. Stability

### Users-22 — No console errors while browsing the Users View (P1)
**Goal:** opening and interacting with the Users View produces no errors.
**Type:** regression
**Preconditions:** error capture on.
**Steps:**
1. Open the Users View. → List renders.
2. Scroll the list; single-click 2–3 users to populate the Context Panel; switch view modes. → Each interaction renders without errors.

**Final check:** no `pageerror`, no `console.error`, no error balloon, no stack in the UI.
**Postconditions/cleanup:** none.
**Selectors:** `watchErrors` / `expectNoErrors`.
