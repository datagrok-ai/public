---
feature: roles
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: []
realizes: []
realized_as:
  - roles.test.ts
related_bugs: []
---

# Roles management — manual test cases

**Instance:** https://dev.datagrok.ai/
**Object under test:** Roles View (Browse > Platform > Roles) — listing, New Role dialog, search,
view modes, context menu, role properties, permission panes, and role assignment (Assigned to →
MANAGE, with "Can assign").
**Type:** UI / manual (staged for a Playwright port)
**Sources:** wiki `govern/access-control/access-control.md` (Global Permissions), `users-and-groups.md`
(roles, "Can assign"); live recon on dev (see `user_groups_tests_context.md`).

---

## Conventions

- **ID:** `Roles-<NN>` — Test Track compatible.
- **Priority:** P1 (critical path) … P3 (edge / cosmetic).
- **Type:** `smoke` / `functional` / `regression` / `negative`.
- Each case: **Goal → Type → Preconditions → Steps (action → expected UI reaction) → Final check →
  Postconditions/cleanup → Selectors**.

**Global preconditions (unless a case says otherwise):**
- Logged in to https://dev.datagrok.ai/ as the **admin** account (`DATAGROK_LOGIN`); app loaded.
- Console-error capture is on for "no errors" cases.
- Browse panel open; Context Panel open (F4).

**What counts as an error:** uncaught `pageerror`; `console.error`; a Datagrok error balloon; a
Dart/JS stack surfaced in the UI or Context Panel.

**Test data:** roles can be created and deleted freely. Created roles use a unique name
`qa_autotest_<timestamp>` and are deleted in cleanup. A role can be assigned to existing subjects
(`All users`, `Administrators`, `Developers`) that are always present on a fresh instance.

---

## 1. Roles View — open and layout

### Roles-01 — Open Roles View from the Browse tree (P1)
**Goal:** Roles View opens from Browse > Platform > Roles.
**Type:** smoke
**Preconditions:** global.
**Steps:**
1. In the Browse tree, expand **Platform**. → Subnodes appear including **Roles**.
2. Click **Roles**. → The Roles View opens; the URL becomes `…/roles`.
3. Inspect the view. → A list of roles is shown (e.g. Administrators, Developers + custom roles); the header has **New Role...**, a search field, view-mode toggles, and a `shown / total` count.

**Final check:** active view type is `roles`; URL contains `/roles`; list and count visible; no errors.
**Postconditions/cleanup:** none.
**Selectors:** `treeNodeByPath(['Platform','Roles'])`, `NEW_ROLE_BUTTON`, `GALLERY_SEARCH`, `GALLERY_COUNTS`.

---

### Roles-02 — Roles View toolbar controls are present (P2)
**Goal:** all expected header controls render.
**Type:** smoke
**Preconditions:** Roles View open.
**Steps:**
1. Inspect the header. → Present: **New Role...** button; **Search roles by name or by #tags** input; view toggles **brief / card / grid**; count `N / N`.

**Final check:** every listed control is visible and enabled.
**Postconditions/cleanup:** none.
**Selectors:** `NEW_ROLE_BUTTON`, `GALLERY_SEARCH`, `VIEW_TOGGLE_BRIEF/CARD/GRID`, `GALLERY_COUNTS`.

---

## 2. Create a role

### Roles-03 — "Create New Role" dialog fields and cancel (P1)
**Goal:** the New Role dialog has the right fields and CANCEL is non-destructive.
**Type:** functional
**Preconditions:** Roles View open.
**Steps:**
1. Click **New Role...**. → The **Create New Role** dialog opens.
2. Inspect fields. → Inputs: **Name**, **Description**; buttons **OK**, **CANCEL**.
3. Click **CANCEL**. → Dialog closes; no role is created; count unchanged.

**Final check:** dialog title `Create New Role`; Name and Description present; count unchanged after CANCEL.
**Postconditions/cleanup:** none.
**Selectors:** `NEW_ROLE_BUTTON`, `.d4-dialog`, `[name="input-host-Name"|"input-host-Description"]`, `.ui-btn`.

---

### Roles-04 — Create a new role end-to-end (P1)
**Goal:** a role created via the dialog appears in the list.
**Type:** functional
**Preconditions:** Roles View open; note total `N`.
**Steps:**
1. Click **New Role...**. → Dialog opens.
2. Enter Name `qa_autotest_<ts>`, Description `created by autotest`. → Fields accept input; OK enabled.
3. Click **OK**. → Dialog closes; the new role is created.
4. Search the Roles View for `qa_autotest_<ts>`. → Exactly one matching role is shown.

**Final check:** the new role is findable by name; total count increased by 1.
**Postconditions/cleanup:** delete the created role (Roles-15).
**Selectors:** `NEW_ROLE_BUTTON`, dialog inputs, `GALLERY_SEARCH`, `GALLERY_COUNTS`.

---

### Roles-05 — New role validation: empty name (P2)
**Goal:** the dialog rejects an empty name.
**Type:** negative
**Preconditions:** Roles View open.
**Steps:**
1. Click **New Role...**. → Dialog opens.
2. Leave **Name** empty and click **OK**. → No role is created: OK is blocked or a validation message/error appears; the dialog stays open.

**Final check:** no role created on empty name; the user is informed.
**Postconditions/cleanup:** CANCEL the dialog.
**Selectors:** dialog inputs, `.ui-btn`, `BALLOON_CONTAINER`.

---

## 3. Search and view modes

### Roles-06 — Search roles by name (P2)
**Goal:** typing a query filters the list and reflects in the URL.
**Type:** functional
**Preconditions:** Roles View open; note total `N`.
**Steps:**
1. Type a known role fragment (e.g. `Admin`) in the search field. → The list narrows; the count drops below `N`; the URL gains `?q=Admin`.
2. Clear the search. → The full list returns; count back to `N / N`; `?q=` removed.

**Final check:** filtered count < total while searching; restored on clear; URL tracks the query.
**Postconditions/cleanup:** clear search.
**Selectors:** `GALLERY_SEARCH`, `GALLERY_COUNTS`, URL assertion `/roles?q=…/`.

---

### Roles-07 — Switch list view modes (P3)
**Goal:** brief / card / grid toggles change the layout.
**Type:** functional
**Preconditions:** Roles View open (brief by default).
**Steps:**
1. Click **Switch to card view**. → The list re-renders as cards.
2. Click **Switch to grid view**. → The list re-renders as a data grid.
3. Click **Switch to brief view**. → Back to the brief list.

**Final check:** each toggle changes the layout and marks itself current (`d4-current`); no errors.
**Postconditions/cleanup:** return to brief view.
**Selectors:** `VIEW_TOGGLE_BRIEF/CARD/GRID`.

---

## 4. Role context menu and properties

### Roles-08 — Role context menu items (P2)
**Goal:** right-clicking a role shows the expected actions.
**Type:** functional
**Preconditions:** Roles View open.
**Steps:**
1. Right-click a custom role (e.g. a `qa_autotest` role). → A context menu opens with: **Properties...**, **Delete**.

**Final check:** the listed items are present. (Built-in roles like Administrators may restrict Delete.)
**Postconditions/cleanup:** close the menu (Esc).
**Selectors:** role card, `.d4-menu-popup .d4-menu-item`, `d4-name="Properties...|Delete"`.

---

### Roles-09 — Edit role properties (rename + description) (P2)
**Goal:** the Properties... dialog edits the role's name and description.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` role exists.
**Steps:**
1. Right-click the role → **Properties...**. → The `<role> Properties` dialog opens with **Name**, **Description** and OK/CANCEL.
2. Change Name to `qa_autotest_<ts>_renamed`, edit the Description. → Fields accept input.
3. Click **OK**. → Dialog closes; the rename is saved.
4. Search the Roles View for `_renamed`. → The renamed role is found.

**Final check:** the role's name/description reflect the edit.
**Postconditions/cleanup:** delete the role (Roles-15).
**Selectors:** context menu `d4-name="Properties..."`, `.d4-dialog`, `[name="input-host-Name"|"input-host-Description"]`.

---

### Roles-10 — Role context-panel info panes (P2)
**Goal:** selecting a role populates the Context Panel with role panes.
**Type:** smoke
**Preconditions:** Roles View open; Context Panel open.
**Steps:**
1. Single-click a role. → The Context Panel shows the role; panes appear: **Actions**, **Assigned to**, **Favorites**, **Global Permissions**, **Permissions**, **Sticky meta**.
2. Expand **Assigned to**. → The subjects (users/groups) the role is assigned to are listed; a **MANAGE** button is present.

**Final check:** the expected panes render; Assigned to lists subjects and shows MANAGE; no errors.
**Postconditions/cleanup:** none.
**Selectors:** role card (click), `.grok-prop-panel .d4-accordion-pane[name="pane-Actions|pane-Assigned-to|pane-Favorites|pane-Global-Permissions|pane-Permissions|pane-Sticky-meta"]`, `ROLE_ASSIGNED_MANAGE_BTN`.

---

## 5. Role assignment

### Roles-11 — Assign a role to a group/user (Assigned to → MANAGE) (P1)
**Goal:** the Assigned to → MANAGE editor assigns the role to a subject.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` role exists; a target subject exists (e.g. `Developers` group or a user).
**Steps:**
1. Single-click the role; expand **Assigned to**; click **MANAGE**. → The `<role> members` dialog opens with a **Search by name or email to add...** input, the current assignee list (each with a **Can assign** toggle), and SAVE/CANCEL.
2. Type a subject name in the search-to-add input and pick a match (user or group). → The subject is added to the assignee list in the dialog.
3. Click **SAVE**. → Dialog closes; the assignment is saved.
4. Reopen **Assigned to**. → The added subject is listed.

**Final check:** the subject is now assigned the role (visible in **Assigned to**).
**Postconditions/cleanup:** delete the role (Roles-15), which removes the assignment.
**Selectors:** `ROLE_ASSIGNED_MANAGE_BTN`, `MEMBERS_ADD_INPUT`, assignee rows, `.ui-btn` (SAVE).

---

### Roles-12 — Remove a role assignment (MANAGE) (P2)
**Goal:** the MANAGE editor removes an assignee.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` role with at least one assignee (from Roles-11).
**Steps:**
1. Open **Assigned to > MANAGE**. → The dialog lists current assignees.
2. Remove an assignee (clear it from the list). → The assignee is removed from the dialog list.
3. Click **SAVE**. → Dialog closes; the removal is saved.
4. Reopen **Assigned to**. → The removed subject is no longer listed.

**Final check:** the assignee is removed.
**Postconditions/cleanup:** delete the role.
**Selectors:** `ROLE_ASSIGNED_MANAGE_BTN`, assignee row remove control, `.ui-btn` (SAVE).

---

### Roles-13 — Grant "Can assign" to a subject (P2)
**Goal:** the per-assignee **Can assign** toggle lets a subject grant the role to others.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` role with an assignee (from Roles-11).
**Steps:**
1. Open **Assigned to > MANAGE**. → Each assignee row has a **Can assign** toggle.
2. Turn on **Can assign** for an assignee. → The toggle becomes active.
3. Click **SAVE**. → Dialog closes; the setting is saved.
4. Reopen **MANAGE**. → The assignee's **Can assign** toggle is still on.

**Final check:** the assignee is persisted with "Can assign" enabled.
**Postconditions/cleanup:** delete the role.
**Selectors:** `ROLE_ASSIGNED_MANAGE_BTN`, assignee row Can-assign toggle, `.ui-btn` (SAVE).

---

### Roles-14 — Edit a role's permissions (Global Permissions / Permissions panes) (P2)
**Goal:** a role's permission set can be edited from the Context Panel.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` role exists.
**Steps:**
1. Single-click the role; expand **Global Permissions** (and **Permissions**). → The pane shows the role's permission set (the exact control — checkboxes vs. an edit dialog — is **to be confirmed on first run**).
2. Toggle one permission on. → The permission becomes selected/saved.
3. Reselect the role / reopen the pane. → The permission change persists.

**Final check:** the toggled permission is reflected for the role.
**Postconditions/cleanup:** revert the permission; delete the role.
**Selectors:** `pane-Global-Permissions`, `pane-Permissions` — **confirm the editing control on first run** (panes rendered lazily during recon).

---

## 6. Delete and stability

### Roles-15 — Delete a role (P1)
**Goal:** Delete removes a role from the list.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` role exists; note total `N`.
**Steps:**
1. Right-click the role → **Delete**. → A confirmation is shown (if any), then the role is removed.
2. Confirm the deletion. → The role disappears from the list.
3. Search for the role name. → No matching role is found; total count is `N - 1`.

**Final check:** the role is gone; count decreased by 1.
**Postconditions/cleanup:** none (this *is* the cleanup for created roles).
**Selectors:** context menu `d4-name="Delete"`, confirmation dialog, `GALLERY_SEARCH`, `GALLERY_COUNTS`.

---

### Roles-16 — No console errors while browsing the Roles View (P1)
**Goal:** opening and interacting with the Roles View produces no errors.
**Type:** regression
**Preconditions:** error capture on.
**Steps:**
1. Open the Roles View. → List renders.
2. Single-click 2–3 roles to populate the Context Panel; expand the Assigned to / Permissions panes; switch view modes. → Each interaction renders without errors.

**Final check:** no `pageerror`, no `console.error`, no error balloon, no stack in the UI.
**Postconditions/cleanup:** none.
**Selectors:** `watchErrors` / `expectNoErrors`.
