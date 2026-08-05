---
feature: groups
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: []
realizes: [views.groups]
realized_as:
  - groups.test.ts
related_bugs: []
---

# Groups management — manual test cases

**Instance:** https://dev.datagrok.ai/
**Object under test:** Groups View (Browse > Platform > Groups) — listing, NEW GROUP dialog, search,
view modes, context menu, group properties, member management (add/remove/admin), nesting, delete.
**Type:** UI / manual (staged for a Playwright port)
**Sources:** wiki `govern/access-control/users-and-groups.md` (Groups), `access-control.md` (permissions);
live recon on dev (see `user_groups_tests_context.md`).

---

## Conventions

- **ID:** `Groups-<NN>` — Test Track compatible.
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

**Test data:** dev already has 40+ groups. Created groups use a unique name `qa_autotest_<timestamp>`
and are **deleted** in cleanup. Unlike users, groups **can** be deleted.

---

## 1. Groups View — open and layout

### Groups-01 — Open Groups View from the Browse tree (P1)
**Goal:** Groups View opens from Browse > Platform > Groups.
**Type:** smoke
**Preconditions:** global.
**Steps:**
1. In the Browse tree, expand **Platform**. → Subnodes appear including **Groups**.
2. Click **Groups**. → The Groups View opens; the URL becomes `…/groups`.
3. Inspect the view. → A list of groups is shown; the header has **NEW GROUP...**, a search field, view-mode toggles, and a `shown / total` count.

**Final check:** active view type is `groups`; URL contains `/groups`; list and count visible; no errors.
**Postconditions/cleanup:** none.
**Selectors:** `treeNodeByPath(['Platform','Groups'])`, `NEW_GROUP_BUTTON`, `GALLERY_SEARCH`, `GALLERY_COUNTS`.

---

### Groups-02 — Groups View toolbar controls are present (P2)
**Goal:** all expected header controls render.
**Type:** smoke
**Preconditions:** Groups View open.
**Steps:**
1. Inspect the header. → Present: **NEW GROUP...** button; **Search groups by name or by #tags** input; view toggles **brief / card / grid**; count `N / N`.

**Final check:** every listed control is visible and enabled.
**Postconditions/cleanup:** none.
**Selectors:** `NEW_GROUP_BUTTON`, `GALLERY_SEARCH`, `VIEW_TOGGLE_BRIEF/CARD/GRID`, `GALLERY_COUNTS`.

---

## 2. Create a group

### Groups-03 — "Create New Group" dialog fields and cancel (P1)
**Goal:** the New Group dialog has the right fields and CANCEL is non-destructive.
**Type:** functional
**Preconditions:** Groups View open.
**Steps:**
1. Click **NEW GROUP...**. → The **Create New Group** dialog opens.
2. Inspect fields. → Inputs: **Name**, **Description**; buttons **OK**, **CANCEL**.
3. Click **CANCEL**. → Dialog closes; no group is created; count unchanged.

**Final check:** dialog title `Create New Group`; Name and Description present; count unchanged after CANCEL.
**Postconditions/cleanup:** none.
**Selectors:** `NEW_GROUP_BUTTON`, `.d4-dialog`, `[name="input-host-Name"|"input-host-Description"]`, `.ui-btn`.

---

### Groups-04 — Create a new group end-to-end (P1)
**Goal:** a group created via the dialog appears in the list.
**Type:** functional
**Preconditions:** Groups View open; note total `N`.
**Steps:**
1. Click **NEW GROUP...**. → Dialog opens.
2. Enter Name `qa_autotest_<ts>`, Description `created by autotest`. → Fields accept input; OK enabled.
3. Click **OK**. → Dialog closes; the new group is created (the group view/Context Panel may open for it).
4. Go back to the Groups View and search `qa_autotest_<ts>`. → Exactly one matching group is shown.

**Final check:** the new group is findable by name; total count increased by 1.
**Postconditions/cleanup:** delete the created group (Groups-14).
**Selectors:** `NEW_GROUP_BUTTON`, dialog inputs, `GALLERY_SEARCH`, `GALLERY_COUNTS`.

---

### Groups-05 — New group dialog: no client-side name validation (P3)
**Goal:** document the Create New Group dialog's validation behavior.
**Type:** functional
**Preconditions:** Groups View open.
**Observed on dev (important):** unlike the **New User** dialog, the **Create New Group** dialog does
**NOT** disable **OK** for an empty name (the OK button stays `enabled`) — there is no client-side
validation. So the automated check only verifies the dialog fields render and **CANCEL** is
non-destructive; it does **not** click OK with an empty name (that would create an unnamed group).
**Steps:**
1. Click **NEW GROUP...**. → Dialog opens with **Name** and **Description**; **OK** is enabled even when Name is empty.
2. Click **CANCEL**. → Dialog closes; no group created; count unchanged.

**Final check:** fields present; CANCEL non-destructive.
**Postconditions/cleanup:** none.
**Selectors:** dialog inputs, `.ui-btn`.

---

## 3. Search and view modes

### Groups-06 — Search groups by name (P1)
**Goal:** typing a query filters the list and reflects in the URL.
**Type:** functional
**Preconditions:** Groups View open; note total `N`.
**Steps:**
1. Type a known group fragment (e.g. `QA`) in the search field. → The list narrows; the count drops below `N`; the URL gains `?q=QA`.
2. Clear the search. → The full list returns; count back to `N / N`; `?q=` removed.

**Final check:** filtered count < total while searching; restored on clear; URL tracks the query.
**Postconditions/cleanup:** clear search.
**Selectors:** `GALLERY_SEARCH`, `GALLERY_COUNTS`, URL assertion `/groups?q=…/`.

---

### Groups-07 — Switch list view modes (P2)
**Goal:** brief / card / grid toggles change the layout.
**Type:** functional
**Preconditions:** Groups View open (brief by default).
**Steps:**
1. Click **Switch to card view**. → The list re-renders as cards.
2. Click **Switch to grid view**. → The list re-renders as a data grid.
3. Click **Switch to brief view**. → Back to the brief list.

**Final check:** each toggle changes the layout and marks itself current (`d4-current`); no errors.
**Postconditions/cleanup:** return to brief view.
**Selectors:** `VIEW_TOGGLE_BRIEF/CARD/GRID`.

---

## 4. Group context menu and properties

### Groups-08 — Group context menu items (P1)
**Goal:** right-clicking a group shows the expected actions.
**Type:** functional
**Preconditions:** Groups View open.
**Steps:**
1. Right-click a group (e.g. a `qa_autotest` group, or `Some demo group`). → A context menu opens with: **Properties...**, **Request membership**, **Chat**, **Delete**.

**Final check:** the listed items are present. (The set may extend with management entries depending on the current user's relation to the group.)
**Postconditions/cleanup:** close the menu (Esc).
**Selectors:** group card, `.d4-menu-popup .d4-menu-item`, `d4-name="Properties...|Request membership|Chat|Delete"`.

---

### Groups-09 — Edit group properties (rename + description) (P1)
**Goal:** the Properties... dialog edits the group's name and description.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` group exists.
**Steps:**
1. Right-click the group → **Properties...**. → The `<group> Properties` dialog opens with **Name**, **Description** and OK/CANCEL.
2. Change Name to `qa_autotest_<ts>_renamed`, edit the Description. → Fields accept input.
3. Click **OK**. → Dialog closes; the rename is saved.
4. Search the Groups View for `_renamed`. → The renamed group is found; the old name no longer matches.

**Final check:** the group's name/description reflect the edit.
**Postconditions/cleanup:** delete the group (Groups-14).
**Selectors:** context menu `d4-name="Properties..."`, `.d4-dialog`, `[name="input-host-Name"|"input-host-Description"]`.

---

### Groups-10 — Group context-panel info panes (P2)
**Goal:** selecting a group populates the Context Panel with group panes.
**Type:** smoke
**Preconditions:** Groups View open; Context Panel open.
**Steps:**
1. Single-click a group. → The Context Panel shows the group; panes appear: **Actions**, **Members**, **Favorites**, **Global Permissions**, **Permissions**, **Sticky meta**.
2. Expand **Members**. → The current members (users and/or nested groups) are listed; a **MANAGE** button is present.

**Final check:** the expected panes render; Members lists members and shows MANAGE; no errors.
**Postconditions/cleanup:** none.
**Selectors:** group card (click), `.grok-prop-panel .d4-accordion-pane[name="pane-Actions|pane-Members|pane-Favorites|pane-Global-Permissions|pane-Permissions|pane-Sticky-meta"]`, `MEMBERS_MANAGE_BTN`.

---

## 5. Member management

### Groups-11 — Add a user to a group (MANAGE) (P1)
**Goal:** the Members → MANAGE editor adds a user to the group.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` group exists; a known user to add (e.g. an existing dev user).
**Steps:**
1. Single-click the group; expand **Members**; click **MANAGE**. → The `<group> members` dialog opens with a **Search by name or email to add...** input, the current member list, and SAVE/CANCEL.
2. Type a user name in the search-to-add input and pick a match. → The user is added to the member list in the dialog.
3. Click **SAVE**. → Dialog closes; membership is saved.
4. Reopen **Members** in the Context Panel. → The added user is listed as a member.

**Final check:** the chosen user is now a member of the group.
**Postconditions/cleanup:** delete the group (Groups-14), which removes the membership.
**Selectors:** `MEMBERS_MANAGE_BTN`, `MEMBERS_ADD_INPUT`, member list rows, `.ui-btn` (SAVE).

---

### Groups-12 — Remove a member from a group (MANAGE) (P2)
**Goal:** the MANAGE editor removes a member.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` group with at least one member (from Groups-11).
**Steps:**
1. Open **Members > MANAGE**. → The dialog lists current members.
2. Remove a member (click its remove control / clear it from the list). → The member is removed from the dialog list.
3. Click **SAVE**. → Dialog closes; the removal is saved.
4. Reopen **Members**. → The removed user is no longer listed.

**Final check:** the member count decreases by one; the removed user is gone.
**Postconditions/cleanup:** delete the group.
**Selectors:** `MEMBERS_MANAGE_BTN`, member row remove control, `.ui-btn` (SAVE).

---

### Groups-13 — Assign a member as group admin (P2)
**Goal:** the per-member **Admin** toggle promotes a member to group admin.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` group with a member (from Groups-11).
**Steps:**
1. Open **Members > MANAGE**. → Each member row has an **Admin** toggle.
2. Turn on **Admin** for a member. → The toggle becomes active for that member.
3. Click **SAVE**. → Dialog closes; the admin assignment is saved.
4. Reopen **MANAGE**. → The member's **Admin** toggle is still on.

**Final check:** the member is persisted as a group admin.
**Postconditions/cleanup:** delete the group.
**Selectors:** `MEMBERS_MANAGE_BTN`, member row Admin toggle, `.ui-btn` (SAVE).

---

### Groups-15 — Nest a group inside another group (P2)
**Goal:** adding a group as a member nests it (members inherit the parent's permissions).
**Type:** functional
**Preconditions:** two `qa_autotest` groups exist — parent `P` and child `C`.
**Steps:**
1. Open parent `P` → **Members > MANAGE**. → The members dialog opens.
2. In the search-to-add input, type child `C`'s name and pick the **group** match. → `C` is added to `P`'s member list (as a group).
3. Click **SAVE**. → Dialog closes; the nesting is saved.
4. Reopen `P`'s **Members**. → Group `C` is listed as a member of `P`.

**Final check:** `C` is a member of `P`; per the model, `C`'s members inherit `P`'s permissions.
**Postconditions/cleanup:** delete both groups.
**Selectors:** `MEMBERS_MANAGE_BTN`, `MEMBERS_ADD_INPUT` (group match), member list.

---

## 6. Membership requests, chat, personal groups, favorites

### Groups-16 — Request membership in a group (P2)
**Goal:** a non-member can request to join a group; the request reaches a group admin.
**Type:** functional
**Preconditions:** **secondary non-admin account** (`DATAGROK_SHARING_LOGIN`); a group the secondary account is not a member of.
**Steps:**
1. As the secondary account, open the Groups View and right-click the target group → **Request membership**. → A confirmation/notification indicates the request was sent.
2. (Optional, as the group admin) Open notifications / the group. → A pending membership request is visible and can be approved/declined.

**Final check:** the request is sent without error; (optionally) it shows up for the group admin.
**Postconditions/cleanup:** the admin declines/cancels the request; remove the member if approved.
**Selectors:** context menu `d4-name="Request membership"`, notification/balloon.
**Observed on dev (NEEDS CLARIFICATION):** clicking **Request membership** as the non-admin surfaced a
**server error** (a balloon with **"Report error"**; a failed `POST` in `DelegatingHttpClient`). It is
unclear whether this is a real bug or an artifact of a duplicate request left by a prior run (the
non-admin has no obvious UI to withdraw its own pending request). The Playwright test is therefore
marked `test.fixme` (`groups.sharing.test.ts`) pending clean-state verification / product
clarification.

---

### Groups-17 — Chat with a group (P3)
**Goal:** the Chat action opens a group chat.
**Type:** functional
**Preconditions:** Groups View open.
**Steps:**
1. Right-click a group → **Chat**. → A chat view/panel for the group opens.

**Final check:** the chat opens without error.
**Postconditions/cleanup:** close the chat.
**Selectors:** context menu `d4-name="Chat"`.

---

### Groups-18 — Each user has a personal group (P3, API-only)
**Goal:** verify a user's personal security group exists.
**Type:** functional
**Observed on dev (important):** personal groups **do exist** (every user's `group` is their personal
group, `isPersonal = true`), but the **Groups View gallery does NOT list them** — searching a user's
name in the Groups View returns **0** results. So there is **no UI** to assert this against; it is
verified at the **API level** instead (not automated in the UI Playwright set).
**Steps (API):**
1. Fetch the user with its `group` (`grok.dapi.users...include('group')`). → `group.isPersonal` is true and `group.name` equals the user's name.

**Final check:** the user's personal group exists and is marked personal (API).
**Postconditions/cleanup:** none. **Do not** delete personal groups.

---

### Groups-19 — Add a group to favorites (P3)
**Goal:** a group can be added to favorites and appears under Favorites.
**Type:** functional
**Preconditions:** Groups View open; Context Panel open.
**Steps:**
1. Single-click a group; in the Context Panel, expand **Favorites** (or use the star). → The favorites control is available.
2. Add the group to favorites. → The group is starred.
3. Open **Browse > My stuff > Favorites**. → The group appears in the favorites list.

**Final check:** the group is present in Favorites.
**Postconditions/cleanup:** remove the group from favorites.
**Selectors:** group card, `pane-Favorites` / `CONTEXT_PANEL_STAR`, Favorites list.

---

## 7. Delete and stability

### Groups-14 — Delete a group (P1)
**Goal:** Delete removes a group from the list.
**Type:** functional
**Preconditions:** a `qa_autotest_<ts>` group exists; note total `N`.
**Steps:**
1. Right-click the group → **Delete**. → A confirmation is shown (if any), then the group is removed.
2. Confirm the deletion. → The group disappears from the list.
3. Search for the group name. → No matching group is found; total count is `N - 1`.

**Final check:** the group is gone; count decreased by 1.
**Postconditions/cleanup:** none (this *is* the cleanup for created groups).
**Selectors:** context menu `d4-name="Delete"`, confirmation dialog, `GALLERY_SEARCH`, `GALLERY_COUNTS`.

---

### Groups-20 — No console errors while browsing the Groups View (P1)
**Goal:** opening and interacting with the Groups View produces no errors.
**Type:** regression
**Preconditions:** error capture on.
**Steps:**
1. Open the Groups View. → List renders.
2. Single-click 2–3 groups to populate the Context Panel; expand the Members pane; switch view modes. → Each interaction renders without errors.

**Final check:** no `pageerror`, no `console.error`, no error balloon, no stack in the UI.
**Postconditions/cleanup:** none.
**Selectors:** `watchErrors` / `expectNoErrors`.
