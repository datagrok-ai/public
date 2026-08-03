# Users, Groups & Roles — test context (handoff)

Companion to `users_manual_tests.md`, `groups_manual_tests.md`, and `roles_manual_tests.md`. Holds everything needed to
continue the work in a fresh session: goal, environment, scope decisions, UI recon results, and
the selector map captured live on dev.

---

## 1. Goal and current stage

- Writing **manual test cases for Users & Groups management** (Browse > Platform > Users / Groups),
  UI-first coverage.
- Final target — **Playwright autotests**. The manual cases are structured for an easy port to code.
- Stage: manual cases drafted (`users_manual_tests.md`, `groups_manual_tests.md`). Next — review,
  then port to Playwright (`@smoke` first, then `functional` / `regression`).

## 2. Environment

- **Target instance:** https://dev.datagrok.ai/
- **Primary account:** `DATAGROK_LOGIN` = `opavlenko+playwright@datagrok.ai` — confirmed member of the
  **Administrators** group, so user/group management actions (create user, block, edit members,
  edit roles, edit any group) are available. Used by default in all cases.
- **Secondary account:** `DATAGROK_SHARING_LOGIN` = `opavlenko+pwsharing@datagrok.ai` — non-admin dev
  account. Needed for cross-account cases: *Request membership*, *invite via URL/password*, and
  negative permission checks. Requires its own `storageState`.
- Credentials live in `playwright-tests/.env`.

## 3. Test-case format

- Markdown, prose (no tables for steps).
- Each case: **Goal → Type → Priority → Preconditions → Steps (action + expected UI reaction) →
  Final check → Postconditions/cleanup → Selectors → ref**.
- **ID:** `Users-<NN>` / `Groups-<NN>` — Test Track compatible.
- **Type:** `smoke` / `functional` / `regression` / `negative` — maps to Playwright tags.
- **Priority:** P1 (critical path) … P3 (edge / cosmetic).
- Atomic steps: one action per line with its expected result. On automation: action → `await`,
  expectation → `expect`. Preconditions → `beforeEach`/fixture; Postconditions → `afterEach`/teardown.

## 4. Scope decisions

- **Cover both views fully**: Users View and Groups View — toolbar controls, NEW dialogs, search,
  view modes, context menus, context-panel info panes, and the management flows (create, block,
  edit memberships/roles, manage members, delete group).
- **UI-first.** Every action that a real user performs goes through the UI. API (`grok.dapi.*`) is
  reserved for setup of preconditions, verification reads, and cleanup — never as a shortcut for the
  action under test.
- **Irreversible user creation drives the two-set split (IMPORTANT).** The platform has **no hard
  delete for users** — a created user can only be **blocked**, not removed. So a "create user" test
  permanently adds a row, which is unacceptable on the shared dev instance. Hence two sets (below):
  the dev set never creates users; the CI/CD set does the real creation on a fresh DB.
- **Sensitive / external-effect actions** (Invite a Friend by email, invite via group password URL):
  cover the dialog/validation in the UI, but do **not** actually send invitations from CI.
- **Roles View** (Browse > Platform > Roles) **is included** — create/edit/delete a role, assign it
  to users/groups (Assigned to → MANAGE, with "Can assign"), and the permission panes.

### Two sets / file layout

The suite exists in **two sets**, differing only in **user creation**:

- **Dev set — `playwright-tests/e2e/user groups/` (reddata).** Runs against the shared **dev**
  instance, so it **never creates users**. Files: `users_manual_tests.md` (existing-users variant),
  `groups_manual_tests.md`, `roles_manual_tests.md`, this context doc. Management cases
  (block, edit Groups.../Roles..., favorites) target **existing `opavlenko<NNN>`** test users and
  **revert** every change in cleanup. Create-user cases (Users-05, Users-07) are **deferred** here.
  Group/Role creation **is** exercised on dev because groups/roles can be deleted (self-clean).
  Also holds the **two-account** scenarios (`Groups-16 Request membership`, permission negatives)
  using the secondary non-admin account.
- **CI/CD set — `public/playwright-public/user groups/` (fresh DB).** Identical coverage **plus the
  active create-user / service-user E2E** (Users-05, Users-07). Currently mirrors
  `users_manual_tests.md` (create variant) + context; the rest is assembled at port time.

**Seeding on CI (important):** on a fresh CI DB the `opavlenko<NNN>` target users do **not** exist, so
the CI/CD set's **create-user step seeds users with the *same names* (`opavlenko<NNN>`)** that the
dev set uses as targets. That single seed both exercises the create flow and provides targets for the
management cases — keeping the two sets structurally aligned (same logins, same case IDs).

## 5. UI recon results (dev, captured live)

### Users View (`/users`, view type `users`)

- **Tree path:** Browse > Platform > **Users**.
- **Header / toolbar:**
  - **NEW** combo button → `User...`, `Service User...`, `Invite a Friend...`.
  - **Search** input, placeholder `Search users by name or by #tags`; supports `#tag`; query is
    reflected in the URL as `?q=<text>`.
  - **View toggles:** brief (`icon-grip-lines`, default) / card (`icon-grip-horizontal`) /
    grid (`icon-table`); **Sort** (`icon-sort-alt`); **Toggle filters** (`icon-filter`);
    **Refresh** (`icon-sync`).
  - **Count:** `.grok-items-view-counts`, e.g. `190 / 190` (shown / total).
- **Dialogs:**
  - **Create new user** — fields: Email, Login, First Name, Last Name; buttons OK / CANCEL.
  - **Create new service user** — field: Login; OK / CANCEL.
  - **Invite a Friend** — field: Email; OK / CANCEL.
- **User context menu (right-click):** `Details`, `Chat`, `Block`, `Groups...`, `Roles...`,
  `ID` (copy entity id), `Add to favorites`.
  - `Groups...` → `<user> groups` dialog: checkbox list of groups, SAVE / CANCEL (edits the user's
    group memberships).
  - `Roles...` → `<user> roles` dialog: role list with a `Can assign` option, SAVE / CANCEL.
  - `Block` → blocks login; the action toggles to unblock for a blocked user.
  - `Details` (or double-click) → opens the **user profile view** (view type `user_edit`).
- **User context-panel info panes:** `Personal`, `Roles`, `Member of`, `Projects`, `Activity`,
  `Clicks`, `Privileges`, `Chats`, `Sticky meta`, `User input`.

### Groups View (`/groups`, view type `groups`)

- **Tree path:** Browse > Platform > **Groups**.
- **Header / toolbar:**
  - **NEW GROUP...** button (ribbon).
  - **Search** input, placeholder `Search groups by name or by #tags`; URL `?q=<text>`.
  - **View toggles:** brief / card / grid.
  - **Count:** `.grok-items-view-counts`, e.g. `41 / 41`.
- **Dialogs:**
  - **Create New Group** — fields: Name, Description; OK / CANCEL.
  - **`<Group> Properties`** (context menu `Properties...`) — fields: Name, Description; OK / CANCEL
    (this is the group edit/rename dialog).
- **Group context menu (right-click):** `Properties...`, `Request membership`, `Chat`, `Delete`.
  (The exact set depends on the current user's relation to the group; admins see management entries.)
- **Group context-panel info panes:** `Actions`, `Members`, `Favorites`, `Global Permissions`,
  `Permissions`, `Sticky meta`.
  - **Members pane** has a **MANAGE** button → `<Group> members` dialog: a `Search by name or email
    to add...` input, the member list with a per-member **Admin** toggle, and SAVE / CANCEL. This is
    the unified add/remove-members + assign-admin editor. Members can be **users or other groups**
    (adding a group nests it — its members inherit the parent group's permissions).
- **Personal groups:** every user automatically has a personal group named after them
  (`isPersonal = true`).

### Roles View (`/roles`, view type `roles`)

- **Tree path:** Browse > Platform > **Roles**.
- **Header / toolbar:** **New Role...** button; search `Search roles by name or by #tags` (URL `?q=`);
  view toggles brief / card / grid; count `N / N`.
- **Dialogs:**
  - **Create New Role** — fields: Name, Description; OK / CANCEL.
  - **`<Role> Properties`** (context menu `Properties...`) — fields: Name, Description; OK / CANCEL.
- **Role context menu (right-click):** `Properties...`, `Delete` (built-in roles may restrict Delete).
- **Role context-panel info panes:** `Actions`, `Assigned to`, `Favorites`, `Global Permissions`,
  `Permissions`, `Sticky meta`.
  - **Assigned to pane** has a **MANAGE** button → `<Role> members` dialog: a `Search by name or
    email to add...` input, the assignee list (users/groups) with a per-assignee **Can assign**
    toggle, and SAVE / CANCEL. This is the role-assignment editor.
  - **Global Permissions / Permissions panes** hold the role's permission set; the exact editing
    control rendered lazily during recon — **confirm on first run**.
- A role can also be assigned from the **Users** view (`Roles...` context menu) and reviewed in a
  user's **Roles** info pane.

## 6. Selector map (to confirm on first Playwright run)

Most controls already have stable selectors from the Browse suite (`browse/selectors.ts`):
tree (`treeNodeByPath(['Platform','Users'])`), context menu (`.d4-menu-popup .d4-menu-item`,
`d4-name="..."`), context panel + info panes (`.grok-prop-panel`, `.d4-accordion-pane[name="pane-X"]`),
dialogs (`.d4-dialog`, `.d4-dialog-title`, `.ui-btn`), inputs (`[name="input-host-<Caption>"]`).

Users/Groups-specific anchors captured live:

| Logical name              | Selector                                                       |
| ------------------------- | -------------------------------------------------------------- |
| `USERS_VIEW`              | view type `users`; URL `/users`                                |
| `GROUPS_VIEW`             | view type `groups`; URL `/groups`                              |
| `NEW_USER_BUTTON`         | `button[name="button-New"]` (combo)                            |
| `NEW_GROUP_BUTTON`        | ribbon button with text `NEW GROUP...`                         |
| `GALLERY_SEARCH`          | `input[placeholder^="Search users"]` / `^="Search groups"`     |
| `GALLERY_COUNTS`          | `.grok-items-view-counts` (`shown / total`)                    |
| `VIEW_TOGGLE_BRIEF`       | `[name="icon-grip-lines"]` (aria `Switch to brief view`)       |
| `VIEW_TOGGLE_CARD`        | `[name="icon-grip-horizontal"]` (aria `Switch to card view`)   |
| `VIEW_TOGGLE_GRID`        | `[name="icon-table"]` (aria `Switch to grid view`)             |
| `SORT_LIST`               | `[name="icon-sort-alt"]` (aria `Sort list`)                    |
| `TOGGLE_FILTERS`          | `[name="icon-filter"]` (aria `Toggle filters`)                 |
| `MEMBERS_MANAGE_BTN`      | text `MANAGE` inside `.grok-prop-panel [name="pane-Members"]`  |
| `MEMBERS_ADD_INPUT`       | `.d4-dialog input[placeholder^="Search by name or email"]`     |
| `NEW_ROLE_BUTTON`         | ribbon button with text `New Role...`                          |
| `ROLE_ASSIGNED_MANAGE_BTN`| text `MANAGE` inside `.grok-prop-panel [name="pane-Assigned-to"]` |
| dialog field by caption   | `[name="input-host-Email" \| "input-host-Login" \| "input-host-First-Name" \| "input-host-Last-Name" \| "input-host-Name" \| "input-host-Description"]` |

## 7. Automation notes

- **baseURL** = `https://dev.datagrok.ai/`; auth via `storageState` from global setup.
- Reuse Browse helpers: `goHome`, `ensureBrowsePanelOpen`, `ensureContextPanelOpen`, `expandTreeGroup`,
  `watchErrors` / `expectNoErrors`, `treeNodeByPath`.
- **Unique names** for created entities: `qa_autotest_<timestamp>`; teardown removes by prefix.
- **Cleanup:** groups — delete via UI then verify gone; users — **block** (no hard delete) and unblock
  in setup if reused. Prefer creating fresh per run over reusing.
- **Wait on UI, not `sleep`:** wait for the new row / dialog / pane, count change, or URL change.
- Cross-account cases (`Request membership`, invite-by-URL) need the secondary non-admin
  `storageState` — model after `*.sharing.test.ts`.

## 8. Open questions

- Should full **create-user** E2E run in CI given users cannot be hard-deleted? (Pollutes the
  instance.) Manual-only, or CI with block-on-cleanup and periodic purge?
- **Invite a Friend** and **invite via group password URL** — verify in UI only, or wire a mailbox /
  second account to assert the invitation actually lands?
- Exact contents of the `ID`/`Copy` user context-menu sub-item (copies id vs. submenu) — confirm on
  first run.
- Whether to cover the **Roles View** (Browse > Platform > Roles) in this set or separately.

## 9. Playwright implementation status (dev set)

Implemented and green on dev (`playwright-tests/e2e/user groups/`, 18 tests, 2 workers):

- `selectors.ts`, `helpers.ts` — gallery/dialog/membership-editor selectors and helpers (reuse Browse
  helpers). Tests are consolidated: one Playwright test covers several manual IDs (title lists them).
- `users.test.ts` (8) — Users-01/02/11, 03/04/06/08, 09, 14, 15/17, 21, 18/19, 20.
- `groups.test.ts` (6) — Groups-01/02/07, 03/04/05/09/14, 06, 08/10, 11/12/13/15.
- `roles.test.ts` (5) — Roles-01/02/07, 03/04/09/15, 06, 08/10, 11/12/13.
- `groups.sharing.test.ts` — Groups-16 (non-admin, `chromium-sharing` project) — **`test.fixme`** (see below).

Key findings baked into the tests/docs:

- **Membership/assignment editor** (`MembershipEditor`, `core/client/d4/.../membership_editor.dart`):
  typing in `.d4-user-selector-input` runs an async search; matches render as `.membership-add-row`
  with `button[name="button-Add"]`; existing rows are `.membership-row` with `button[name="button-Remove"]`
  and a per-row Admin / "Can assign" checkbox. Backs user Groups.../Roles..., group Members→MANAGE,
  role Assigned to→MANAGE. Real keyboard input (`pressSequentially`) is required to trigger the search.
- **No UI unblock** — blocking is `Block` → "Are you sure?" → YES (async `status=blocked`); the gallery
  menu does NOT flip to Unblock and the profile has no control. Unblock only via `POST /api/users/unblock`
  (used for cleanup). Documented UI-first exception.
- **Create New Group / Role** do NOT disable OK on empty name (no client-side validation) — unlike New User.
- **Personal groups** exist (API) but are NOT listed in the Groups View → Groups-18 is API-only.
- **Request membership** (non-admin) surfaced a **server error** (failed POST, "Report error" balloon) —
  `groups.sharing.test.ts` is `test.fixme` pending clean-state verification / product clarification.
- Tests that mutate the same `opavlenko<NNN>` user run `test.describe.configure({ mode: 'serial' })`.

Run (**1 worker** — these tests mutate shared `opavlenko45` and create/delete entities, so they must
not run in parallel): `npx playwright test users.test.ts groups.test.ts roles.test.ts --project=chromium --workers=1`
(sharing: `--project=chromium-sharing`). `REUSE_AUTH=1` reuses `e2e/.auth.json`. Full run ≈ **3.7 min**, 18 tests green.

Stability notes:
- **Freshly-created entities lag the server search index** — a new group/role card briefly appears
  from cache then a server search returns 0 and removes it. Always locate just-created entities with
  `searchAndWaitCard()` (search + Refresh + retry until the card sticks), not a bare `searchGallery` +
  `toBeVisible`. This was the main flake source.
- View toggles scoped to `.grok-items-view-toggle`; clear search via the X reset icon (not `fill('')`).
- A worker-scoped reused page (`closeAll` between tests) was tried for speed but proved fragile with
  retries (shared page breaks on a timeout → cascade) — reverted to isolated per-test pages.
