---
feature: browse
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: []
realizes: [views.browse]
realized_as:
  - apps.test.ts
  - apps_matrix.test.ts
  - ctxpanel.test.ts
  - dash.test.ts
  - db.test.ts
  - demo_apps.test.ts
  - fav.test.ts
  - files.test.ts
  - filter.test.ts
  - modelhub.test.ts
  - mystuff.test.ts
  - nav.test.ts
  - permission.sharing.test.ts
  - platform.test.ts
  - route.test.ts
  - search.test.ts
  - section18.test.ts
  - spaces.test.ts
  - splitmodes.test.ts
  - tree.test.ts
  - view.test.ts
related_bugs: [GROK-16261, GROK-16857, GROK-17664, GROK-17896, GROK-17922, GROK-19628, GROK-19638, GROK-19688, GROK-19689, GROK-19690, GROK-19691, GROK-19692, GROK-19740, GROK-19802, GROK-19844, GROK-19847, GROK-19848, GROK-19934, GROK-19965, GROK-20032]
---

# Datagrok Browse — manual test cases

**Instance:** https://dev.datagrok.ai/
**Object under test:** Browse (navigation, tree, view modes), all tree sections
**Type:** UI / manual (draft to be ported to Playwright)
**Sources:** docs `datagrok.ai/help/datagrok/navigation/`; Jira GROK (Test Track, Browse bugs)

---

## Conventions

- Case **ID**: `Browse-<Area>-<NN>` — format compatible with Test Track (cf. `Scripts-Edit-6`).
- **Priority:** P1 (critical path) … P3 (edge / cosmetic).
- **Type:** `smoke` (basic operability), `functional` (behavior verification), `regression` (tied to a specific bug), `negative` (invalid/empty input, lack of permissions).
- Each case: **Goal** → **Type** → **Preconditions** → **Steps** (action + expected UI reaction) → **Final check** → **Postconditions/cleanup**.
- Steps are atomic: one action per line with an expected result. When automating, action → `await` call, expectation → `expect`. Preconditions → `beforeEach`/fixture; "Postconditions" → `afterEach`/teardown.
- Criteria are measurable: where possible, specific values are given (timings, sets of menu items, exact state) rather than "correct/relevant".
- Selectors are not yet fixed (we will capture them on the first run on dev, see "Notes for automation"). In the text, elements are named by their UI label/role.
- `ref:` — related bug/ticket in Jira.

**Global preconditions (unless a case states otherwise):**
- The user is logged in at https://dev.datagrok.ai/ under an account with rights to the demo data; the start page is open.
- Console-error interception is enabled (see below).
- **Context Panel is open** in all Browse tests (toggle F4). If the panel closes during a test — reopen it. For section 18 ("tree nodes without errors"), the Context Panel must be **fully expanded** (Expand all info panes). This requirement does **not** apply to the future matrix sets for Demo / Apps.

**Test accounts (see `playwright-tests/.env`):**
- `DATAGROK_LOGIN` = `opavlenko+playwright@datagrok.ai` — primary, **admin** (used by default in all cases).
- `DATAGROK_SHARING_LOGIN` = `opavlenko+pwsharing@datagrok.ai` — second dev account, **without admin rights**. Used in cases involving sharing and permission-restriction checks (Tree-07, Platform-03, MyStuff-04, etc.). For such cases we prepare a separate `storageState`.

**What we consider an error (for the "no errors" checks):**
- `pageerror` (unhandled page exception);
- `console.error` in the browser console;
- a datagrok balloon with error text (red/error notification);
- a Dart/JS stack surfaced in the UI or Context Panel.
Any one of these is enough to consider the case failed.

**Test data:** dev has demo files (Files > Demo), demo DBs (Databases), and demo dashboards (Dashboards). The cases rely on them. For entities created during a test, use a unique name prefix (e.g. `qa_autotest_<timestamp>`) to simplify cleanup.

---

## 1. Browse — general navigation and controls

### Browse-Nav-01 — Open the Browse panel from the Sidebar (P1)
**Goal:** Browse opens and contains the mandatory tree sections.
**Type:** smoke
**Preconditions:** global.
**Steps:**
1. On the Sidebar (leftmost strip), click the **Browse** icon. → The Browse panel with the tree expands on the left.
2. Inspect the top level of the tree. → The **mandatory** nodes are present: **My stuff**, **Spaces**, **Apps**, **Files**, **Dashboards**, **Databases**, **Platform**.

**Final check:** the Browse panel is visible; all 7 mandatory sections are present at the top level (additional ones may exist — we do not treat them as an error); no errors.
**Selectors:** `SIDEBAR_BROWSE_ICON`, `BROWSE_PANEL_TAB_PANE`, `BROWSE_HEADER`, `treeGroupByName('My stuff' | 'Spaces' | 'Apps' | 'Files' | 'Databases' | 'Platform')`, `treeItemByName('Dashboards')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Nav-02 — Hide the Browse panel by clicking again (P2)
**Goal:** toggle behavior of the Browse icon.
**Type:** functional
**Preconditions:** the Browse panel is open.
**Steps:**
1. On the Sidebar, click the **Browse** icon again. → The Browse panel hides, the central workspace expands.
2. Click the **Browse** icon once more. → The panel reappears with the same tree.

**Final check:** after the first click the panel is not visible; after the second — it is visible again.
**Selectors:** `SIDEBAR_BROWSE_ICON`, `BROWSE_PANEL_TAB_PANE`, `SIDEBAR_TAB_SELECTED` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Nav-03 — The Home button opens the Home Page (P1)
**Goal:** navigate to the Home Page from the Browse Top Menu.
**Type:** smoke
**Preconditions:** the Browse panel is open.
**Steps:**
1. In the Browse panel's Top Menu, click the **Home** icon. → The **Home Page** view opens in the center.
2. Inspect the Home Page. → A global search bar is present at the top; below it — the widgets area.

**Final check:** the active view is the Home Page; the search field is present and available for input.
**Selectors:** `BROWSE_HEADER_HOME`, `HOME_VIEW_HANDLE`, `HOME_GLOBAL_SEARCH_INPUT` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Nav-04 — Collapse all collapses the tree (P2)
**Goal:** bulk collapse of nodes.
**Type:** functional
**Preconditions:** ≥2 nodes are expanded in the tree at different levels.
**Steps:**
1. In the Browse Top Menu, click **Collapse all**. → All nodes collapse down to the top level.

**Final check:** no child node is expanded; only the 7 top-level sections are visible; the count of expanded nodes = 0.
**Selectors:** `BROWSE_HEADER_COLLAPSE_ALL`, `TREE_EXPAND_ARROW_EXPANDED` (to check for absence) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Nav-05 — Locate current object highlights the node (P2)
**Goal:** navigate to the node of the currently open object.
**Type:** functional
**Preconditions:** a demo file is open in the center; its tree node is collapsed/not highlighted.
**Steps:**
1. In the Browse Top Menu, click **Locate current object**. → The tree expands to the node of the current object and highlights it.

**Final check:** the current object's node is visible within the scroll area and marked as selected; the node name matches the name of the open object.
**Selectors:** `BROWSE_HEADER_LOCATE`, `treeNodeByName(currentObjectName)`, `treeNodeSelected` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Nav-06 — Refresh tree pulls in server-side changes (P1)
**Goal:** refresh the tree without reloading the page, without losing state.
**Type:** functional
**Preconditions:** the tree is open; there is a way to create an object on the server outside the current tab (recommended via API/second context in setup).
**Steps:**
1. Expand the target section (e.g. My Files) and record its contents + which nodes are expanded. → Initial state captured.
2. Create a new object on the server in that section (via API/second context). → The object is not yet visible in the tree.
3. In the Browse Top Menu, click **Refresh tree**. → The tree is re-read.

**Final check:** the new object appears in the tree without a page reload; the page did not reload (same URL and session); no errors.
**Selectors:** `BROWSE_HEADER_REFRESH`, `treeNodeByName(createdObjectName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** delete the created object.

---

### Browse-Nav-07 — Import file opens the data in a Table View (P2)
**Goal:** import a local file from the Top Menu.
**Type:** functional
**Preconditions:** a local CSV with known columns/row count is at hand.
**Steps:**
1. In the Browse Top Menu, click **Import file**. → The system file-picker dialog opens.
2. Select the local CSV. → The file loads.

**Final check:** a Table View opens; the number of rows/columns matches the source CSV; no errors.
**Selectors:** `BROWSE_HEADER_IMPORT_FILE`, `viewTabHandle('Table')`, `STATUS_BAR_VIEW_PANEL` (Rows/Columns) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the Table View; delete the temporary table if needed.

---

### Browse-Nav-08 — Import text (P3)
**Goal:** open the text-import view.
**Type:** smoke
**Preconditions:** global.
**Steps:**
1. In the Browse Top Menu, click **Import text**. → The **Import text** view opens with an input field.

**Final check:** the Import text view is active; the input field is available.
**Selectors:** `BROWSE_HEADER_IMPORT_TEXT`, `viewTabHandle('Import text')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view.

---

### Browse-Nav-09 — All Browse header icons work (P1, matrix)
**Goal:** no icon in the `.panel-titlebar.disable-selection.grok-browse-header` header leads to an error.
**Type:** smoke + matrix
**Preconditions:** the Browse panel is open.
**Steps:**
1. Enumerate all visible clickable icons in `.grok-browse-header` (record their order and labels/title). → A complete list is obtained (fix the exact set on the first UI-recon run; expected minimum: Home, Import file, Import text, Refresh tree, Collapse all, Locate current object).
2. For each icon in turn: click → verify the expected reaction (opening the corresponding view / triggering the action) → return to the initial state (if needed, e.g. close the Import text view). → At each step: no console.error, no pageerror, no error balloon.

**Final check:** the list of actually visible icons ⊇ the expected minimum; for each icon the click leads to a correct reaction; none causes an error.
**Selectors:** `BROWSE_HEADER`, `BROWSE_HEADER_ICONS`, `BROWSE_HEADER_HOME`, `BROWSE_HEADER_IMPORT_FILE`, `BROWSE_HEADER_IMPORT_TEXT`, `BROWSE_HEADER_REFRESH`, `BROWSE_HEADER_COLLAPSE_ALL`, `BROWSE_HEADER_LOCATE`, `BROWSE_HEADER_CLOSE_PANEL` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close all views opened during the test.

> Cases Nav-03..08 cover the **specific behavior** of each icon. Nav-09 is a matrix smoke — "all icons click without errors", useful when new icons are added.

---

## 2. Browse tree — working with nodes

### Browse-Tree-01 — Expand and collapse a node (P1)
**Goal:** expand / collapse a node by clicking the arrow.
**Type:** smoke
**Preconditions:** the tree is open, the **Files** section is collapsed.
**Steps:**
1. Click the expand arrow to the left of the **Files** node. → The node expands, child items appear.
2. Inspect the child nodes. → My files, Spaces, App Data, Demo are visible; the arrow is in the "expanded" state.
3. Click the collapse arrow of the **Files** node. → The node collapses, child items are hidden.

**Final check:** after expanding, the child nodes are visible and the arrow is "expanded"; after collapsing, the child nodes are absent from the DOM/invisible and the arrow is "collapsed".
**Selectors:** `treeGroupByName('Files')`, `treeExpandArrow('Files')`, `TREE_EXPAND_ARROW_EXPANDED`, `treeNodeChildren(page, 'Files')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Tree-02 — Keyboard navigation in the tree (P2)
**Goal:** control the tree from the keyboard — navigation, expand/collapse, activation.
**Type:** functional
**Preconditions:** focus is set on a tree node (click a node with children, children collapsed).
**Steps:**
1. Press **↓**. → Focus moves to the next visible node.
2. Press **↑**. → Focus returns to the original node.
3. Press **→**. → The node expands, child items appear.
4. Press **←**. → The node collapses.
5. On an item node (item, not group), press **Enter** (or **Space**). → The item opens in the center (like a single click).

**Final check:** ↑↓ switch the selection (at any moment exactly one node is selected); → expands, ← collapses; Enter/Space opens the item.
**Selectors:** `treeNodeContainer(page, name)`, `treeNodeSelected`, `TREE_EXPAND_ARROW_EXPANDED`, keys `ArrowUp` / `ArrowDown` / `ArrowLeft` / `ArrowRight` / `Enter` / `Space` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Tree-03 — The tree remembers expand/collapse state (P1)
**Goal:** preserve the tree state when switching views.
**Type:** regression
**Preconditions:** specific nodes are expanded (e.g. Files > Demo and Databases), the rest collapsed.
**Steps:**
1. Record the list of expanded nodes. → Initial state captured.
2. Open another view (Home Page). → Browse loses focus.
3. Return to Browse (the Sidebar icon). → The tree is displayed again.

**Final check:** the set of expanded nodes is identical to the one captured in step 1; nothing extra was expanded or collapsed. ref: GROK-16261
**Selectors:** `TREE_EXPAND_ARROW_EXPANDED` (to count the expanded nodes), `SIDEBAR_BROWSE_ICON`, `HOME_VIEW_HANDLE` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** collapse the tree (Collapse all) to keep the next test clean.

---

### Browse-Tree-04 — The side menu does not "stick" expanded (P2)
**Goal:** correct state of nested nodes when collapsing/opening the side menu.
**Type:** regression
**Preconditions:** several nested nodes are expanded (record which).
**Steps:**
1. Collapse the vertical side menu. → The menu hides.
2. Open the side menu again. → The menu is displayed.

**Final check:** the set of expanded nodes matches what it was before collapsing; no forced expansion of all nested items occurred. ref: GROK-19802
**Selectors:** `SIDEBAR`, `TREE_EXPAND_ARROW_EXPANDED` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** Collapse all.

---

### Browse-Tree-05 — Context menu on right-click (P1)
**Goal:** availability and composition of the context menu on an item.
**Type:** functional
**Preconditions:** a demo dashboard (entity) is visible in the tree.
**Steps:**
1. Right-click the dashboard. → The context menu opens.
2. Inspect the items. → A minimal set is present: **Open**, **Add to Favorites**, **Share**, **Rename**, **Delete** (verify and fix the exact list on the first run).

**Final check:** the menu is open; the expected items are present and clickable; for a file (non-entity) the set differs — no Share/Rename as for an entity (cross-check with Browse-Fav-06).
**Selectors:** `CONTEXT_MENU`, `CONTEXT_MENU_ITEM`, `contextMenuItem(page, 'Open' | CONTEXT_MENU_ADD_FAVORITES | CONTEXT_MENU_SHARE | CONTEXT_MENU_RENAME | CONTEXT_MENU_DELETE)`, `treeNodeByName(dashboardName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the menu (Esc).

---

### Browse-Tree-06 — Reorganize the tree with drag-and-drop (P2)
**Goal:** move an item in the tree by dragging.
**Type:** functional
**Preconditions:** there is a movable test entity and a valid target position; record the initial location.
**Steps:**
1. Grab the item and start dragging. → Available drop zones are highlighted.
2. Release the item in the new position. → The item moves.

**Final check:** the item is in the new position after a tree refresh; it is absent from the initial position; no errors.
**Selectors:** `treeNodeByName(source)`, `treeNodeByName(target)`, `TREE_NODE_DROP`, Playwright API: `dragTo` or `mouse.down/move/up` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** return the item to its initial position (or recreate the state via a fixture).

---

### Browse-Tree-07 — A node without permissions opens without an error (P2, negative)
**Goal:** clicking an object without permissions does not break the UI.
**Type:** negative
**Preconditions:** logged in as `DATAGROK_SHARING_LOGIN` (non-admin); the tree contains an object whose content this user cannot access (prepare via a fixture/second user).
**Steps:**
1. Click the object without permissions. → A preview stub / a "insufficient permissions" message opens.

**Final check:** a meaningful access-restriction message is shown (not a blank screen and not an error stack); no console errors. ref: GROK-17922
**Selectors:** `treeNodeByName(noAccessObjectName)`, `BALLOON_CONTAINER`, `ERROR_BALLOON` — see `selectors.md` / `selectors.ts`. **Account:** `DATAGROK_SHARING_LOGIN`.
**Postconditions/cleanup:** none required.

---

## 3. Browsing mode vs Persistent view

### Browse-View-01 — Single (unpinned) vs double (persistent) click (P1)
**Goal:** browsing mode — a single click creates an unpinned view (gets replaced), a double click creates a persistent one (not replaced).
**Type:** functional
**Preconditions:** the tree contains three demo files A, B, C.
**Steps:**
1. Single-click A. → A opens as an unpinned view.
2. Single-click B. → View A is replaced by view B; the number of open tabs did not grow.
3. Double-click C. → C opens in a persistent view; B stays open.
4. Single-click A. → A opens, C stays open (not displaced).

**Final check:** unpinned is replaced on the next single click; persistent is preserved after clicking another item; at the end A (current) and C (pinned) are open.
**Selectors:** `treeNodeByName(A | B | C)`, `viewTabHandle(viewName)`, `VIEW_TAB`, `VIEW_TAB_SELECTED`, `VIEW_TAB_CLOSE` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close all views.

---

### Browse-View-02 — Editing makes a view persistent (P2)
**Goal:** unpinned → persistent on any edit.
**Type:** functional
**Preconditions:** an unpinned view is open (single click).
**Steps:**
1. Make a change in the view (e.g. a change in the table). → The view is pinned as persistent.
2. Single-click another tree item. → The edited view is not replaced.

**Final check:** the edited view stays open after clicking another item.
**Selectors:** `viewTabHandle(viewName)`, `treeNodeByName(otherEntity)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view without saving.

---

### Browse-View-03 — Pin pins a view manually (P2)
**Goal:** pinning via pin in tabbed mode.
**Type:** functional
**Preconditions:** tabbed mode is on; an unpinned view is open.
**Steps:**
1. Hover over the blue left border of the tab. → The pin control appears.
2. Click **pin**. → The view is pinned.
3. Single-click another tree item. → The pinned view is not replaced.

**Final check:** the view is marked as pinned and stayed after clicking another item.
**Selectors:** `VIEW_TAB_PIN` (⏳ TODO), `viewTabHandle(viewName)`, `treeNodeByName(otherEntity)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view.

---

### Browse-View-04 — A new browsing session after persistent (P2)
**Goal:** return to browsing mode without converting a persistent view back.
**Type:** functional
**Preconditions:** a persistent view is open.
**Steps:**
1. On the Sidebar, click the **Browse** icon. → The Browse panel is active.
2. Single-click an item. → The item opens in a new unpinned view.

**Final check:** the previously pinned view stays persistent; the new item is in browsing mode (replaced on the next click).
**Selectors:** `SIDEBAR_BROWSE_ICON`, `viewTabHandle(viewName)`, `treeNodeByName(entity)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the open views.

---

### Browse-View-05 — Grouping same-type views with a badge (P3)
**Goal:** a counter for same-type views on the Sidebar.
**Type:** functional
**Preconditions:** 2 views of the same type are open (two tables).
**Steps:**
1. Look at the Sidebar. → Views of the same type are grouped under one icon with a badge.

**Final check:** the badge shows a number = 2; when one view is closed, the badge becomes = 1.
**Selectors:** `SIDEBAR_VIEW_BADGE` (⏳ TODO), `viewTabHandle(viewName)`, `VIEW_TAB_CLOSE` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close all open views.

---

## 4. My stuff

### Browse-MyStuff-01 — Composition of the My stuff section (P1)
**Goal:** presence of subgroups.
**Type:** smoke
**Preconditions:** the tree is open.
**Steps:**
1. Expand **My stuff**. → Subgroups appear.
2. Inspect the subgroups. → The **mandatory** ones are present: Recent, Favorites, Shared with me, My files.
3. Additionally: sub-nodes for specific categories (`My dashboards`, `My scripts`) are absent on the current dev — content is flat under My stuff. If they return as virtual subgroups in the future — add them to the list.

**Final check:** the listed mandatory subgroups are present; user content (dashboards, scripts) is either nested in its own subgroup or flat under My stuff — both variants are valid.
**Selectors:** `treeGroupByName('My stuff')`, `treeExpandArrow('My stuff')`, `treeNodeByPath(['My-stuff', 'Recent'|'Favorites'|'Shared-with-me'|'My-files'])` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-MyStuff-02 — Recent contains recent objects (P2)
**Goal:** freshness of Recent.
**Type:** functional
**Preconditions:** the target demo file has not yet been opened in the current session.
**Steps:**
1. Open a specific demo file (record its name). → The file is open.
2. Expand **My stuff > Recent**. → The list updates.

**Final check:** this exact file is present in Recent (by name), preferably at the top of the list.
**Selectors:** `treeNodeByName('Recent')`, `treeNodeChildren(page, 'Recent')`, `treeNodeByName(fileName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-MyStuff-03 — Favorites shows the favorites (P1)
**Goal:** display of favorite entities.
**Type:** functional
**Preconditions:** exactly one known entity is in favorites (create it in setup).
**Steps:**
1. Open **My stuff > Favorites**. → The favorites list is shown.

**Final check:** the prepared entity is present in the list by name; the number of items matches the prepared one.
**Selectors:** `treeNodeByName('Favorites')`, `treeNodeChildren(page, 'Favorites')`, `treeNodeByName(entityName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** remove the entity from favorites.

---

### Browse-MyStuff-04 — Shared with me is grouped by user (P2)
**Goal:** grouping of shared content.
**Type:** functional
**Preconditions:** in the fixture, user `DATAGROK_SHARING_LOGIN` (U) shares a known object O with user `DATAGROK_LOGIN`; the test runs as `DATAGROK_LOGIN`.
**Steps:**
1. Open **My stuff > Shared with me**. → The shared objects are shown.
2. Inspect the grouping. → There is a group with a header = the name/email of user U, and object O inside it.

**Final check:** object O is located in a group whose header corresponds to U.
**Selectors:** `treeNodeByName('Shared with me')`, `treeNodeChildren(page, 'Shared with me')`, `treeGroupByName(userU)`, `treeNodeByName(objectO)` — see `selectors.md` / `selectors.ts`. **Fixture accounts:** `DATAGROK_SHARING_LOGIN` (shares), `DATAGROK_LOGIN` (test).
**Postconditions/cleanup:** revoke the share if needed.

---

### Browse-MyStuff-05 — Add to Favorites from My Files (P1)
**Goal:** add to favorites from the context menu in My Files.
**Type:** regression
**Preconditions:** My Files contains a known file/entity that is not in favorites.
**Steps:**
1. Right-click an item in **My Files**. → The context menu opens.
2. Check the **Add to Favorites** item. → The item is present and enabled.
3. Choose **Add to Favorites**. → The item is added to favorites.

**Final check:** the item appears in My stuff > Favorites (by name). ref: GROK-19848
**Selectors:** `treeNodeByName('My files')`, `treeNodeByName(itemName)`, `CONTEXT_MENU`, `contextMenuItem(page, CONTEXT_MENU_ADD_FAVORITES)`, `treeNodeByName('Favorites')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** remove the item from favorites.

---

### Browse-MyStuff-06 — An object without a project lands in My stuff (P2)
**Goal:** default placement of a new object.
**Type:** functional
**Preconditions:** global.
**Steps:**
1. Create/upload an object with a unique name, without choosing a project/space. → The object is saved.
2. Open **My stuff**. → The list updates.

**Final check:** the object with this name is found in My stuff; it is not present in other projects/Spaces.
**Selectors:** `treeGroupByName('My stuff')`, `treeNodeByName(createdObjectName)`, `treeGroupByName('Spaces')` (verify the object is absent) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** delete the created object.

---

## 5. Favorites

### Browse-Fav-01 — Add to favorites from the context menu (P1)
**Goal:** Add to Favorites from the tree.
**Type:** functional
**Preconditions:** the tree contains a known entity that is not in favorites.
**Steps:**
1. Right-click the entity. → The context menu opens.
2. Choose **Add to Favorites**. → The entity is added.

**Final check:** the entity is present in Favorites (by name).
**Selectors:** `treeNodeByName(entityName)`, `CONTEXT_MENU`, `contextMenuItem(page, CONTEXT_MENU_ADD_FAVORITES)`, `treeNodeByName('Favorites')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** remove from favorites.

---

### Browse-Fav-02 — Toggle favorite via the star in the Context Panel (P1)
**Goal:** a roundtrip add → remove via the star icon in the Context Panel.
**Type:** functional
**Preconditions:** an entity is open, the Context Panel is visible; the entity is **not** in favorites (the star is neutral).
**Steps:**
1. In the Context Panel, next to the name, click the neutral star. → The star turns orange.
2. Open **My stuff > Favorites**. → The entity is present in the list.
3. Return to the entity's Context Panel and click the orange star. → The star turns neutral.
4. Open **My stuff > Favorites** again. → The entity is absent from the list.

**Final check:** the roundtrip works: after the first click — the star is orange and the entity is in Favorites; after the second click — the star is neutral and the entity is absent from Favorites.
**Selectors:** `CONTEXT_PANEL_STAR`, `treeNodeByName('Favorites')`, `treeNodeChildren(page, 'Favorites')`, `treeNodeByName(entityName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required (the test rolls back its own state).

---

### Browse-Fav-03 — Open the Favorites Panel from the Sidebar (P2)
**Goal:** access favorites from the Sidebar.
**Type:** smoke
**Preconditions:** ≥1 known entity is in favorites.
**Steps:**
1. On the Sidebar, click the **Favorites** icon. → The Favorites Panel opens.

**Final check:** the panel shows a list containing the prepared entity.
**Selectors:** `SIDEBAR_FAVORITES_ICON`, `SIDEBAR_FAVORITES_ICON_BY_DATA`, `treeNodeByName(entityName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** remove from favorites if needed.

---

### Browse-Fav-04 — Add to favorites by dragging (P3)
**Goal:** favorite via drag-and-drop.
**Type:** functional
**Preconditions:** the Favorites Panel/view is open; the tree contains a known entity that is not in favorites.
**Steps:**
1. Drag the entity from the tree into the Favorites Panel/view. → The drop zone is highlighted.
2. Release. → The entity is added.

**Final check:** the entity appears in favorites (by name).
**Selectors:** `treeNodeByName(entityName)`, `SIDEBAR_FAVORITES_ICON`, `TREE_NODE_DROP`, Playwright `dragTo` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** remove from favorites.

---

### Browse-Fav-05 — A file/cell value cannot be added to favorites (P3, negative)
**Goal:** favorites restricted to entities only.
**Type:** negative
**Preconditions:** a standalone file (non-entity) and a table with values are open.
**Steps:**
1. Right-click the standalone file. → The context menu has no Add to Favorites item (or it is disabled).
2. In the file's Context Panel. → The star icon is absent/unavailable.
3. Right-click/context on a cell value. → There is no way to add it to favorites.

**Final check:** neither a file nor a cell value can be added to favorites; no errors.
**Selectors:** `CONTEXT_MENU`, `contextMenuItem(page, CONTEXT_MENU_ADD_FAVORITES)` (verify absence), `CONTEXT_PANEL_STAR` (verify absence) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

## 6. Apps

> **Detailed coverage (opening each application, regression matrix)** is moved into a separate autotest set (`apps_matrix.test.ts`). For each application in the list, **3 checks** are run: (1) single click → the view opens without errors; (2) double click → the view opens without errors; (3) if the view has a `Run` / `OK` / `Submit` button — click it with default parameters, wait for the load, verify no errors. If there is no button — step (3) is skipped.
>
> This section covers only the basic smoke checks of the Apps section in the Browse tree.

### Browse-Apps-01 — The Apps list loads without errors and in a reasonable time (P1)
**Goal:** display Apps + absence of abnormal delays.
**Type:** smoke + regression
**Preconditions:** ≥1 application is installed.
**Steps:**
1. Record the time t0 and expand **Apps**. → The list starts loading.
2. Wait for the list to render fully (t1). → The list is displayed.
3. Inspect the list. → The list is non-empty; the expected applications are present (e.g. Tutorials — fix the exact set on the first run).

**Final check:** the list is non-empty; (t1 − t0) < 5 seconds (fine-tune the threshold with the team; it used to be ~40 s); no errors. ref: GROK-20032
**Selectors:** `treeGroupByName('Apps')`, `treeExpandArrow('Apps')`, `treeNodeChildren(page, 'Apps')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Apps-02 — Open an application without errors (P1)
**Goal:** launch an application from the tree — smoke on one reference application.
**Type:** smoke
**Preconditions:** the Apps list contains a reference application (choose a stable one, e.g. Tutorials).
**Steps:**
1. Click the reference application. → The application view opens.

**Final check:** the application view is open (the title/content matches the application); no errors.
**Selectors:** `treeNodeByName('Tutorials')` (or another reference), `viewTabHandle(appName)`, `VIEW_CLOSE` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the application view.

> Coverage of **every application** in the installed list is in `apps_matrix.test.ts` (a parameterized test over the application list).

---

### Browse-Apps-03 — Application tooltip/details without errors (P2)
**Goal:** a correct hover tooltip.
**Type:** regression
**Preconditions:** Apps contains an application (e.g. the Chem area).
**Steps:**
1. Hover over the application item. → A tooltip appears.
2. Inspect the tooltip and the Context Panel details. → No error messages, no "null"/stack.

**Final check:** the tooltip and details contain meaningful text; no errors. ref: GROK-19638
**Selectors:** `treeNodeByName(appName)`, Playwright `hover`, `CONTEXT_PANEL`, `CONTEXT_PANEL_INNER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

## 7. Files

> **Detailed coverage (opening each demo file from Files > Demo)** is moved into a separate autotest set (`demo_files_matrix.test.ts`). This section has only the basic smoke / regression cases.

### Browse-Files-01 — Composition of the Files section (P1)
**Goal:** presence of the file-browser sections.
**Type:** smoke
**Preconditions:** the tree is open.
**Steps:**
1. Expand **Files**. → The sections appear.
2. Inspect the sections. → My files, App Data are present (on dev — **two nodes** with the same name App Data, both are considered OK), Demo.

**Final check:** all expected Files sections are present (My files, App Data ≥1, Demo). There is currently no "Spaces" section inside Files on dev — this is normal.
**Selectors:** `treeGroupByName('Files')`, `treeNodeChildren(page, 'Files')`, `treeNodeByName('My files' | 'App Data' | 'Demo')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Files-02 — Preview of a tabular demo file (P1)
**Goal:** open the preview of a tabular file.
**Type:** smoke
**Preconditions:** Files > Demo contains a known tabular file with a known row/column count.
**Steps:**
1. Open **Files > Demo**. → The contents are shown.
2. Click the tabular file. → An interactive preview opens (a temporary dataframe).

**Final check:** the preview displays the data; the row/column count matches the expected; no errors.
**Selectors:** `treeNodeByName('Demo')`, `treeNodeByName(tableFileName)`, `viewTabHandle(...)`, `STATUS_BAR_VIEW_PANEL` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the preview.

---

### Browse-Files-03 — Open a folder as a folder view (P2)
**Goal:** navigate through folders.
**Type:** functional
**Preconditions:** Files contains a folder with a known item count.
**Steps:**
1. Click a folder in the Files tree. → A folder view opens.

**Final check:** the folder view shows the contents; the item count matches the expected.
**Selectors:** `treeNodeByName(folderName)`, `viewTabHandle(folderName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Files-04 — A new file in S3/Spaces after Refresh (P1)
**Goal:** refresh the tree for external storage.
**Type:** regression
**Preconditions:** an S3/Spaces storage is connected; there is a way to add a file on the server (via API/second context).
**Steps:**
1. Record the contents of the target folder. → Initial state captured.
2. Add a new file on the server outside the tab. → The file is not yet visible.
3. Click **Refresh tree**. → The tree is re-read.

**Final check:** the new file (by name) appeared in the tree. ref: GROK-19844
**Selectors:** `BROWSE_HEADER_REFRESH`, `treeNodeByName(newFileName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** delete the added file.

---

### Browse-Files-05 — A shared folder appears under My stuff > Shared with me (P2)
**Goal:** verify that a folder shared by a user correctly lands in `My stuff > Shared with me`, grouped by the sharer's name.
**Type:** regression
**Preconditions:** another user has shared a folder with the current user (fixture: dev already has several — `Olesia Pavlenko`, `Andrew Golovko`, ...).
**Steps:**
1. Expand **My stuff > Shared with me**. → Sub-nodes appear — the names of the users who shared something.
2. Expand the first user node (e.g. `Olesia Pavlenko`). → Inside, the specific shared entities/folders are visible.

**Final check:** under Shared with me there are ≥1 user groups; the first user group has ≥1 child item. ref: GROK-19847.

**Note on the visual icon:** on the current dev the folder icon is the same for shared and regular folders (`grok-icon grok-ds-files grok-ds`). The sharing indication is structural (position in the tree), not graphical. If the icon changes — update the test.

**Selectors:** `treeGroupByName('My stuff')`, `treeGroupByName('Shared with me')`, `[name^="tree-My-stuff---Shared-with-me---<user>"]` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Files-06 — Download a file from the context menu (P2)
**Goal:** download via the context menu.
**Type:** functional
**Preconditions:** Files contains a known file.
**Steps:**
1. Right-click the file. → The context menu opens.
2. Choose the download item. → The download starts.

**Final check:** the file is downloaded; the name/extension matches the source; size > 0.
**Selectors:** `treeNodeByName(fileName)`, `CONTEXT_MENU`, `contextMenuItem(page, CONTEXT_MENU_DOWNLOAD)`, Playwright `page.on('download', ...)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** delete the downloaded file from downloads (in the autotest — from the temporary download folder).

---

## 8. Dashboards

### Browse-Dash-01 — List of available dashboards (P1)
**Goal:** display dashboards according to permissions.
**Type:** smoke
**Preconditions:** ≥1 accessible demo dashboard exists.
**Steps:**
1. Click the **Dashboards** node in the tree. → The Dashboards node is an **item list** (does not expand); the click opens the **Dashboards** view with a list.

**Final check:** the Dashboards view is open; the list contains the known dashboards `chemical_space_demo`, `demo-datagrok-api`; those not accessible by permissions are not shown.
**Selectors:** `treeItemByName('Dashboards')`, `viewTabHandle('Dashboards')`, `LIST_SEARCH_INPUT` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view.

---

### Browse-Dash-02 — Open a dashboard (P1)
**Goal:** open a dashboard in a Table View.
**Type:** smoke
**Preconditions:** the list contains a known dashboard.
**Steps:**
1. Click the dashboard. → A Table View opens.

**Final check:** the Table View displays the data; the title matches the dashboard; ≥1 viewer is visible; no errors.
**Selectors:** `treeNodeByName(dashboardName)`, `viewTabHandle(dashboardName)`, `STATUS_BAR_VIEW_PANEL` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view.

---

### Browse-Dash-03 — The Context Panel matches the selected dashboard (P3)
**Goal:** absence of duplicated content in the Context Panel.
**Type:** regression
**Preconditions:** there are 2 different dashboards A and B with distinguishable content.
**Steps:**
1. Click A. → The Context Panel shows A's data.
2. Click B. → The Context Panel updates to B's data.

**Final check:** the Context Panel content for A and B differs (e.g. different chat/description text); no identical text. ref: GROK-19934
**Selectors:** `treeNodeByName(A | B)`, `CONTEXT_PANEL`, `CONTEXT_PANEL_INNER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

## 9. Databases

### Browse-DB-01 — List of connected DBs (P1)
**Goal:** display the top-level databases.
**Type:** smoke
**Preconditions:** a connected demo DB exists (known name).
**Steps:**
1. Expand **Databases**. → The list of databases is shown.

**Final check:** the known demo DB is present at the top level.
**Selectors:** `treeGroupByName('Databases')`, `treeNodeByName('Postgres')` (or another known DB) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-DB-02 — Expand a DB down to tables and columns (P1)
**Goal:** the connection → tables → columns → queries hierarchy.
**Type:** functional
**Preconditions:** there is a demo DB with a known table and columns.
**Steps:**
1. Expand the database. → Connections appear.
2. Expand down to tables. → The known table is visible.
3. Expand the table. → Its columns are visible; saved queries are present (if any).

**Final check:** the hierarchy expands down to columns; the table/column names match the expected.
**Selectors:** `treeNodeByName(dbProvider | connectionName | tableName | columnName)`, `treeExpandArrow(name)`, `TREE_EXPAND_ARROW_EXPANDED` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** Collapse all.

---

### Browse-DB-03 — Schema Browser from the context menu (P2)
**Goal:** launch the Schema Browser.
**Type:** functional
**Preconditions:** there is a connection in Databases.
**Steps:**
1. Right-click the connection. → The context menu opens.
2. Choose **Browse schema**. → The Schema Browser opens.

**Final check:** the Schema Browser is open; tables and relations are displayed; no errors.
**Selectors:** `treeNodeByName(connectionName)`, `CONTEXT_MENU`, `contextMenuItem(page, CONTEXT_MENU_BROWSE_SCHEMA)`, `viewTabHandle('Schema Browser')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the Schema Browser.

---

### Browse-DB-04 — Open a saved query (P1)
**Goal:** open the Query Editor.
**Type:** functional
**Preconditions:** a saved query exists.
**Steps:**
1. Click the query. → The Query Editor opens.
2. Inspect the Context Panel. → The SQL code and parameters are shown.

**Final check:** the Query Editor is open; the Context Panel shows the SQL (non-empty) and a parameters block.
**Selectors:** `treeNodeByName(queryName)`, `viewTabHandle(queryName)`, `CONTEXT_PANEL`, `infoPaneByName(page, 'SQL')` (exact pane name — ⏳ to be confirmed) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the Query Editor.

---

### Browse-DB-05 — Info panes for a table (P2)
**Goal:** context info panes for a table.
**Type:** functional
**Preconditions:** there is a demo DB with a table.
**Steps:**
1. Click the table. → The Context Panel updates.
2. Inspect the info panes. → The panes are present: metadata (Details), a dynamic content preview, a Run/query action.

**Final check:** the listed info panes are present; the preview loads the table rows; no errors.
**Selectors:** `treeNodeByName(tableName)`, `CONTEXT_PANEL`, `CONTEXT_PANEL_ACCORDION_PANE`, `infoPaneByName(page, 'Details')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-DB-06 — Postgres > CHEMBL > Browse > Summary without errors (P3)
**Goal:** regression of a specific path.
**Type:** regression
**Preconditions:** a Postgres connection with CHEMBL is available.
**Steps:**
1. Expand Postgres > CHEMBL > Browse. → The items are visible.
2. Click the last item **Summary**. → A view opens.

**Final check:** Summary opens; no errors. ref: GROK-16857
**Selectors:** `treeNodeByName('Postgres')`, `treeNodeByName('CHEMBL')`, `treeNodeByName('Browse')`, `treeNodeByName('Summary')`, `viewTabHandle('Summary')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

## 10. Platform

### Browse-Platform-01 — Composition of the Platform section (P1)
**Goal:** presence of the instance-configuration sections.
**Type:** smoke
**Preconditions:** admin/power user rights.
**Steps:**
1. Expand **Platform**. → The sections appear.
2. Inspect the sections. → The plugins, functions, and access-control sections are present (fix the exact list on the first run).

**Final check:** the expected configuration sections are present and clickable.
**Selectors:** `treeGroupByName('Platform')`, `treeNodeChildren(page, 'Platform')`, `treeNodeByName('Plugins' | 'Functions' | 'Users' | 'Groups' | 'Roles' | 'Admin')` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** Collapse all.

---

### Browse-Platform-02 — Open a Platform section (P2)
**Goal:** open the gallery/view of a section.
**Type:** functional
**Preconditions:** Platform is expanded.
**Steps:**
1. Open the section with plugins or functions. → The gallery/view opens.

**Final check:** the section view is open (the title matches); no errors.
**Selectors:** `treeNodeByName(sectionName)`, `viewTabHandle(sectionName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view.

---

### Browse-Platform-03 — Platform restriction by permissions (P3, negative)
**Goal:** hiding/restricting Platform for a regular user.
**Type:** negative
**Preconditions:** logged in as `DATAGROK_SHARING_LOGIN` (non-admin) via a separate `storageState`.
**Steps:**
1. Open Browse as the non-admin user. → The tree is displayed.
2. Check the visibility/availability of **Platform**. → The section is hidden or its internal actions are blocked.

**Final check:** unavailable sections/actions do not execute; on an access attempt — a correct restriction without an error stack.
**Selectors:** `treeGroupByName('Platform')` (verify absence/disabled), `BALLOON_CONTAINER` (no error balloon) — see `selectors.md` / `selectors.ts`. **Account:** `DATAGROK_SHARING_LOGIN`.
**Postconditions/cleanup:** switch back to the primary account (admin).

---

## 11. Filtering & Search within Browse

### Browse-Filter-01 — The filter panel without empty properties (P1)
**Goal:** relevance of the filter property list.
**Type:** regression
**Preconditions:** an entity gallery is open (Files or Apps).
**Steps:**
1. Open the filter panel. → The property list is shown.
2. Inspect the properties. → There are no properties with an empty name/without values.

**Final check:** the panel has no empty/nameless properties; all listed properties are applicable to the current entity type. ref: GROK-19691
**Selectors:** `FILTER_TOGGLE`, `FILTER_PANEL` (⏳ TODO), `filterProperty(page, name)` (⏳ TODO) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** reset the filters.

---

### Browse-Filter-02 — The CreatedRecently filter for Files (P1)
**Goal:** filter files by creation date + reset.
**Type:** regression
**Preconditions:** the Files list is open; there are recent and old files.
**Steps:**
1. Apply the **CreatedRecently** filter. → The list is rebuilt.
2. Remove the filter. → The list is restored.

**Final check:** with the filter on, only recently created files are shown; after removal — the full list as before the filter. ref: GROK-19689
**Selectors:** `filterProperty(page, 'CreatedRecently')`, `LIST_SEARCH_INPUT` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** ensure the filter is removed.

---

### Browse-Filter-03 — The UsedByMe filter for Apps/Dockers (P1)
**Goal:** filter by usage by the current user.
**Type:** regression
**Preconditions:** the Apps or Dockers list is open; some items have been used by the user, some have not.
**Steps:**
1. Apply the **UsedByMe** filter. → The list is rebuilt.
2. Remove the filter. → The list is restored.

**Final check:** with the filter on, only used items are shown; after removal — the full list. ref: GROK-19688
**Selectors:** `filterProperty(page, 'UsedByMe')`, `LIST_SEARCH_INPUT` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** ensure the filter is removed.

---

### Browse-Filter-04 — The Users filter by group without errors (P2)
**Goal:** filter users by group.
**Type:** regression
**Preconditions:** the Users list is open; there is a group with a known number of members.
**Steps:**
1. Apply the filter by group. → The list is rebuilt.

**Final check:** only the group members are shown (the count matches the expected); no balloons/errors. ref: GROK-19692
**Selectors:** `filterProperty(page, 'Group')`, `LIST_SEARCH_INPUT`, `BALLOON_CONTAINER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** reset the filter.

---

### Browse-Filter-05 — Filter and search working together (P1)
**Goal:** filter + search simultaneously.
**Type:** regression
**Preconditions:** a list is open where ≥1 item is known to match both the filter AND the search string.
**Steps:**
1. Apply the filter. → The list narrows by the filter.
2. Enter the search string. → The list narrows further.

**Final check:** the result = the intersection of filter and search (the expected known item is present, non-matching ones are absent). ref: GROK-19690
**Selectors:** `filterProperty(page, name)`, `LIST_SEARCH_INPUT` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** reset the filter and search.

---

### Browse-Filter-06 — Empty filter/search result (P2, negative)
**Goal:** a correct empty state.
**Type:** negative
**Preconditions:** a list with search capability is open.
**Steps:**
1. Enter a deliberately non-existent string in the search (e.g. `zzz_no_match_zzz`). → The list is rebuilt.

**Final check:** an empty result / a "nothing found" message is shown; no errors and no "stuck" previous list.
**Selectors:** `LIST_SEARCH_INPUT`, `LIST_EMPTY_STATE` (⏳ TODO) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** clear the search field.

---

## 12. Routing (URL state)

### Browse-Route-01 — Restore a dashboard by URL (P2)
**Goal:** open a view via a saved link.
**Type:** functional
**Preconditions:** a dashboard is open; access is granted by permissions.
**Steps:**
1. Copy the URL from the address bar. → The link is obtained.
2. Open the link in a new tab. → Datagrok loads the view.

**Final check:** the same dashboard opens (the title matches); the state (e.g. the selected Table View) is restored.
**Selectors:** Playwright `page.goto(url)`, `viewTabHandle(dashboardName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the tab/view.

---

### Browse-Route-02 — Open a file/folder by URL (P2)
**Goal:** direct access to a file/folder.
**Type:** functional
**Preconditions:** a URL of the form `/files/{namespace}.{share}/{path}` to an existing nested folder is known.
**Steps:**
1. Open the nested folder's URL. → Datagrok navigates to the path.

**Final check:** exactly the specified folder is open (at the correct depth); its contents are displayed.
**Selectors:** Playwright `page.goto('/files/{ns}.{share}/{path}')`, `viewTabHandle(folderName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-Route-03 — Parameterized query by URL (P3)
**Goal:** run a query from a link with parameters.
**Type:** functional
**Preconditions:** there is a saved parameterized query with a known expected result.
**Steps:**
1. Open the URL `/q/{ns}.{conn}.{query}?param=value`. → The query runs with the parameter.

**Final check:** the result matching the parameter value is shown (e.g. the row count/values match the expected).
**Selectors:** Playwright `page.goto('/q/{ns}.{conn}.{query}?param=value')`, `viewTabHandle(queryName)`, `STATUS_BAR_VIEW_PANEL` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view.

---

### Browse-Route-04 — An invalid URL does not break the application (P3, negative)
**Goal:** resilience to a broken link.
**Type:** negative
**Preconditions:** global.
**Steps:**
1. Open a deliberately non-existent entity URL (e.g. `/p/no.such_project/none`). → Datagrok handles the request.

**Final check:** a "not found" message / redirect to a safe screen is shown; the application does not crash, there is no error stack in the UI.
**Selectors:** Playwright `page.goto(brokenUrl)`, `BALLOON_CONTAINER`, `APP_LOADED` (the application did not crash) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

## 13. Spaces — smoke

> Full Spaces coverage (rename, moving, permission inheritance, data loss) is in a separate Spaces test set. Here — only a basic check of the node in Browse.

### Browse-Spaces-01 — Click Spaces + display of available spaces (P1)
**Goal:** the Spaces node clicks without errors and shows the spaces available to the user.
**Type:** smoke
**Preconditions:** the Browse tree is open; ≥1 available space is known.
**Steps:**
1. Click the **Spaces** node (expand via the arrow). → The node expands; the available spaces are shown.

**Final check:** the known space is present; those not accessible by permissions are not shown; no errors (neither console.error, nor pageerror, nor a balloon).
**Selectors:** `treeGroupByName('Spaces')`, `treeExpandArrow('Spaces')`, `treeNodeChildren(page, 'Spaces')`, `treeNodeByName(spaceName)`, `BALLOON_CONTAINER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

## 14. Model Hub / Compute in Browse — smoke

> Deep testing of Compute/Model Hub is out of the Browse scope. Here — a smoke over the points where errors historically occurred (GROK-19965, 19740, 19628, 17896, 17770).

### Browse-ModelHub-01 — Open the Model Catalog from Apps (P1)
**Goal:** launch the Model Catalog without errors.
**Type:** regression
**Preconditions:** Compute / Model Hub is installed.
**Steps:**
1. Open **Apps > Compute > Model Catalog** (Model Hub). → The model-catalog view opens.

**Final check:** the view is open; the model list is displayed; no errors. ref: GROK-17896, GROK-17664
**Selectors:** `treeNodeByName('Apps')`, `treeNodeByName('Compute')`, `treeNodeByName('Model Catalog')` (or 'Model Hub'), `viewTabHandle(...)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view.

---

### Browse-ModelHub-02 — Single click on a model in the tree (P1)
**Goal:** model preview without errors on click/hover.
**Type:** regression
**Preconditions:** Model Hub is expanded in the tree, models are visible.
**Steps:**
1. Single-click a model. → The model opens in a preview.
2. Hover over the model in the tree. → A tooltip is shown.

**Final check:** the preview is rendered; the tooltip contains meaningful text; no errors. ref: GROK-19740
**Selectors:** `treeNodeByName(modelName)`, Playwright `hover`, `CONTEXT_PANEL_INNER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the preview.

---

### Browse-ModelHub-03 — Double click vs right-click → Run (P2)
**Goal:** regression of the double-click bug + control of the working alternative.
**Type:** regression
**Preconditions:** Model Hub is expanded in the tree.
**Steps:**
1. Double-click a model. → The model opens in a persistent view.
2. (Control) Right-click another model → **Run**. → The model launches.

**Final check:** the double click gives no errors (it historically crashed); right-click → Run still works. ref: GROK-19965
**Selectors:** `treeNodeByName(modelName)`, Playwright `dblclick`, `CONTEXT_MENU`, `contextMenuItem(page, CONTEXT_MENU_RUN)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the open views.

---

### Browse-ModelHub-04 — Expand Uncategorized models (P2)
**Goal:** regression of the group-expansion bug.
**Type:** regression
**Preconditions:** the tree contains an **Uncategorized models** group.
**Steps:**
1. Expand the **Uncategorized** group. → The group expands, the models are shown.

**Final check:** the group is expanded, the models are visible; no errors. ref: GROK-19628
**Selectors:** `treeGroupByName('Uncategorized models')`, `treeExpandArrow('Uncategorized models')`, `TREE_EXPAND_ARROW_EXPANDED` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** collapse the group.

---

## 15. Global search (Home Page)

> Global search is available from the bar on the Home Page. The documentation describes several mechanisms: search by entity metadata, by tags, by identifiers, NLP. Here — basic coverage of metadata/tags.

### Browse-Search-01 — Find an entity by name and navigate to it (P1)
**Goal:** find an entity by name and open it from the results.
**Type:** functional
**Preconditions:** on the Home Page; a known entity with a unique name is known (demo dashboard/file).
**Steps:**
1. Enter the name of the known entity in the search bar. → The widgets area is replaced by the results.
2. Inspect the results. → The searched-for entity is present.
3. Click the found entity. → The corresponding view opens.

**Final check:** the entity is present in the results; after the click the correct view is open, the title/content matches the selected result; no errors.
**Selectors:** `HOME_VIEW_HANDLE`, `HOME_GLOBAL_SEARCH_INPUT`, `HOME_SEARCH_RESULTS` (⏳ TODO), `viewTabHandle(entityName)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the view; clear the search bar.

---

### Browse-Search-02 — Search by tag (P2)
**Goal:** filter by tag.
**Type:** functional
**Preconditions:** on the Home Page; a tag assigned to ≥1 entity is known (e.g. `#demo`).
**Steps:**
1. Enter `#demo` in the search bar. → Results are shown.

**Final check:** all results are marked with the `#demo` tag; irrelevant ones are not shown.
**Selectors:** `HOME_GLOBAL_SEARCH_INPUT`, `HOME_SEARCH_RESULTS` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** clear the search bar.

---

### Browse-Search-03 — Empty search result (P2, negative)
**Goal:** a correct empty state of global search.
**Type:** negative
**Preconditions:** on the Home Page.
**Steps:**
1. Enter a deliberately non-existent string (`zzz_no_match_zzz`). → Search results are shown.

**Final check:** "nothing found" / an empty results block is shown; no errors; the widgets do not "break".
**Selectors:** `HOME_GLOBAL_SEARCH_INPUT`, `HOME_SEARCH_RESULTS`, `HOME_SEARCH_EMPTY_STATE` (⏳ TODO), `HOME_WIDGETS` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** clear the search bar (the widgets return).

---

## 16. Context Panel controls

> The Context Panel is the right-hand auxiliary screen of Browse. Toggle via F4; in the header — Back/Forward, Clone and detach, Collapse all / Expand all, Favorites.

### Browse-CtxPanel-01 — Toggle the Context Panel via F4 (P2)
**Goal:** show/hide the panel with a hotkey.
**Type:** functional
**Preconditions:** the Context Panel is visible.
**Steps:**
1. Press **F4**. → The Context Panel hides.
2. Press **F4** again. → The Context Panel is displayed again.

**Final check:** F4 toggles the Context Panel visibility; after pressing again the content is the same.
**Selectors:** the `F4` key, `CONTEXT_PANEL` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** leave the panel visible.

---

### Browse-CtxPanel-02 — Context Panel updates when the object changes (P1)
**Goal:** the panel content follows the current object.
**Type:** functional
**Preconditions:** the Context Panel is visible; the tree has two distinguishable objects A and B.
**Steps:**
1. Click A. → The panel shows A's data (name/Details correspond to A).
2. Click B. → The panel updates to B's data.

**Final check:** the panel title/Details correspond to the current object; the content differs for A and B.
**Selectors:** `treeNodeByName(A | B)`, `CONTEXT_PANEL`, `CONTEXT_PANEL_INNER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-CtxPanel-03 — Back/Forward through viewed objects (P2)
**Goal:** navigation through the Context Panel history.
**Type:** functional
**Preconditions:** objects A and then B were viewed in sequence (the panel shows B).
**Steps:**
1. In the Context Panel header, click **Back**. → The panel returns to A.
2. Click **Forward**. → The panel shows B again.

**Final check:** Back returns to the previous object (A), Forward — to the next one (B).
**Selectors:** `CONTEXT_PANEL_BACK`, `CONTEXT_PANEL_FORWARD`, `treeNodeByName(A | B)`, `CONTEXT_PANEL_INNER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

### Browse-CtxPanel-04 — Collapse all / Expand all info panes (P3)
**Goal:** bulk collapse/expand of info panes.
**Type:** functional
**Preconditions:** an object with ≥2 info panes is open.
**Steps:**
1. In the header, click **Collapse all**. → All info panes collapse.
2. Click **Expand all**. → All info panes expand.

**Final check:** after Collapse all the number of expanded panes = 0; after Expand all all are expanded.
**Selectors:** `CONTEXT_PANEL_COLLAPSE_ALL`, `CONTEXT_PANEL_EXPAND_ALL`, `CONTEXT_PANEL_ACCORDION_PANE` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** none required.

---

## 17. Split view / viewing modes

> A view can be split (Hamburger → Split right/down) and the display modes can be switched from the Status Bar: Tabbed / Simple / Presentation.

### Browse-Split-01 — Split right divides the workspace (P2)
**Goal:** split the central area.
**Type:** functional
**Preconditions:** one view is open (e.g. a demo table).
**Steps:**
1. Click the **Hamburger** in the view's Top Menu. → The menu opens.
2. Choose **Split right**. → The workspace is split into two areas horizontally.
3. Click another object in the tree. → It opens in one of the areas.

**Final check:** two views are visible side-by-side at the same time; both are interactive.
**Selectors:** `VIEW_HAMBURGER` (⏳ TODO), `VIEW_HAMBURGER_SPLIT_RIGHT` (⏳ TODO), `treeNodeByName(otherEntity)`, `viewTabHandle(...)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close one of the views (return to single mode).

---

### Browse-Split-02 — Split down divides the workspace vertically (P3)
**Goal:** vertical split.
**Type:** functional
**Preconditions:** one view is open.
**Steps:**
1. Hamburger → **Split down**. → The workspace is split into top/bottom.

**Final check:** two views are arranged vertically; both are interactive.
**Selectors:** `VIEW_HAMBURGER` (⏳ TODO), `VIEW_HAMBURGER_SPLIT_DOWN` (⏳ TODO), `viewTabHandle(...)` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close one of the views.

---

### Browse-Modes-01 — Switch Tabbed / Simple / Presentation (P2)
**Goal:** change the view display mode.
**Type:** functional
**Preconditions:** ≥2 views are open.
**Steps:**
1. On the Status Bar, switch to **Simple**. → Views are switched via a dropdown selector (tabs are hidden).
2. Switch to **Presentation**. → Only the active view is visible, the rest of the UI is hidden.
3. Switch back to **Tabbed**. → Views are displayed as tabs again.

**Final check:** each mode gives the expected layout; switching does not lose open views and does not cause errors.
**Selectors:** `STATUS_BAR`, `STATUS_BAR_MODE_TABBED` / `STATUS_BAR_MODE_SIMPLE` / `STATUS_BAR_MODE_PRESENTATION` (⏳ TODO — exact icon mapping) — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** return to Tabbed; close the extra views.

---

## 18. Browse — tree nodes open without errors (matrix)

> **Goal of the section:** ensure that **clicking each node of the Browse tree does not lead to an error** (no console.error, no pageerror, no error balloon, no stack in the Context Panel).
>
> **Global preconditions for section 18:**
> - The Context Panel is **open and fully expanded** (Expand all info panes).
> - If the Context Panel closes during a test — reopen and expand it.
> - On every action, check `assertNoErrors(page)` (console / pageerror / balloon).

---

### Browse-Node-MyStuff-01 — My stuff: click each sub-node + the first item (P1, matrix)
**Goal:** all My stuff sub-nodes open without errors; for non-empty ones — the first item in the list opens without errors.
**Type:** smoke + matrix
**Preconditions:** global + Context Panel expanded.
**Steps (for each node in the list):**
- Mandatory nodes: `Recent`, `Favorites`, `Shared with me`, `My dashboards`, `My scripts`, `My files`.
- Optional nodes: `my connections`, `my others`, `my spaces`, `my tables` — if present on dev, perform the same steps; if absent — skip.

For each node:
1. Click the node in the tree. → The node expands / a view with a list opens.
2. If there is a child item in the tree — click the first one. Otherwise, if a view with a list opened — click the first item in the view. Otherwise skip the step. → A preview / view of the selected item opens.
3. On every step, `assertNoErrors(page)`. → Not a single error.

**Final check:** for each present node the steps ran without errors.
**Selectors:** `treeGroupByName('My stuff')`, `treeNodeByName(subNodeName)`, `treeNodeChildren(page, subNodeName)`, `CONTEXT_PANEL`, `CONTEXT_PANEL_EXPAND_ALL`, `BALLOON_CONTAINER`, `ERROR_BALLOON` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close the open views; collapse My stuff (Collapse all).

---

### Browse-Node-Spaces-01 — Spaces: click the parent node + expand (P2)
**Goal:** the Spaces node clicks and expands without errors.
**Type:** smoke
**Preconditions:** global + Context Panel expanded.
**Steps:**
1. Click **Spaces**. → The node expands; the available spaces are visible (on dev: `P`, `PW-Gen-SrcCopy-20`, `PW-Gen-SrcLink-19`, `PW-Gen-Tgt-19`, `PW-Gen-Tgt-20`).
2. `assertNoErrors(page)`. → No errors.

**Final check:** Spaces expanded (or reported "empty" if there is nothing); no errors. Going inside each space is **not needed**.
**Selectors:** `treeGroupByName('Spaces')`, `treeExpandArrow('Spaces')`, `TREE_EXPAND_ARROW_EXPANDED`, `BALLOON_CONTAINER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** collapse Spaces.

---

### Browse-Node-Apps-* — Apps: matrix walkthrough of all applications (P1, matrix → separate file)

> **Coverage of each application** is in `apps_matrix.test.ts` (separately). See the note at the beginning of section 6.

---

### Browse-Node-Files-Demo-* — Files > Demo: matrix walkthrough (P1, matrix → separate file)

> **Coverage of each demo file** is in `demo_files_matrix.test.ts` (separately). See the note at the beginning of section 7.

In section 18 — only a basic check: open Demo, open 1–2 subfolders, open 1–2 files from them (automator's choice). See also Browse-Files-02 (preview of a tabular demo file).

---

### Browse-Node-Dashboards-01 — Dashboards: open the view + specific dashboards (P1)
**Goal:** the Dashboards node opens a view with a list; the known demo dashboards open without errors.
**Type:** smoke
**Preconditions:** global + Context Panel expanded.
**Steps:**
1. Click **Dashboards** in the tree. → The Dashboards view with a list opens.
2. Find and open `chemical_space_demo`. → The corresponding view opens.
3. `assertNoErrors(page)`. → No errors.
4. Return to Dashboards and open `demo-datagrok-api`. → A view opens.
5. `assertNoErrors(page)`. → No errors.

**Final check:** both dashboards are open, without errors.
**Selectors:** `treeItemByName('Dashboards')`, `viewTabHandle('Dashboards')`, `LIST_SEARCH_INPUT`, `treeNodeByName('chemical_space_demo' | 'demo-datagrok-api')`, `BALLOON_CONTAINER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close all open views.

---

### Browse-Node-DB-01 — Databases: click each top-level provider (P1, matrix)
**Goal:** each DB provider in the tree clicks/expands without errors.
**Type:** smoke + matrix
**Preconditions:** global + Context Panel expanded.
**Steps:** Expand **Databases**. For each provider in the list:
- On dev the following are present: Access, Athena, BigQuery, Cassandra, ClickHouse, DB2, Databricks, Denodo, Firebird, HBase, Hive, Hive2, Impala, MLFlow, MS SQL, MariaDB, MongoDB, MySQL, Neo4j, Neptune, ODATA, Oracle, PI, Postgres, Redshift, SAP HANA, Snowflake, Sparql, Teradata, Vertica, Virtuoso.

For each:
1. Click the provider node. → The node expands / shows connections in the Context Panel.
2. `assertNoErrors(page)`. → No errors.

**Final check:** each provider clicked without errors.
**Selectors:** `treeGroupByName('Databases')`, `treeNodeByName(providerName)`, `treeExpandArrow(providerName)`, `CONTEXT_PANEL`, `BALLOON_CONTAINER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** Collapse all.

---

### Browse-Node-DB-Postgres-01 — Postgres: in-depth check (P1)
**Goal:** verify opening connections, tables, queries, and the Schema Browser in Postgres.
**Type:** functional + regression
**Preconditions:** global + Context Panel expanded; Databases > Postgres is expanded.
**Steps:**
1. Expand `Postgres > CHEMBL`. Expand down to tables, click 1–2 tables. Expand one table down to columns, click 1–2 columns. → The hierarchy expands; the Context Panel updates for each; no errors. (Browse-DB-02, DB-05.)
2. Right-click the `CHEMBL` connection → **Browse schema**. → The Schema Browser opens; no errors. (Browse-DB-03.)
3. Close the Schema Browser. Expand `Postgres > CHEMBL > Browse`, click the **Summary** item. → Summary opens; no errors. ref: GROK-16857. (Browse-DB-06.)
4. Expand `Postgres > Datagrok`. Open the saved queries: `World`, `NorthwindTest`. For each: the Query Editor opens; the Context Panel shows non-empty SQL and a parameters block. → No errors. (Browse-DB-04.)
5. (Optional) The same with `Postgres > ChemblTest` — expand down to a couple of tables.

**Final check:** all steps without errors; ref: GROK-16857.
**Selectors:** `treeNodeByName('Postgres' | 'CHEMBL' | 'ChemblTest' | 'Datagrok' | 'Browse' | 'Summary' | 'World' | 'NorthwindTest')`, `CONTEXT_MENU`, `contextMenuItem(page, CONTEXT_MENU_BROWSE_SCHEMA)`, `viewTabHandle('Schema Browser' | 'Summary' | queryName)`, `CONTEXT_PANEL`, `infoPaneByName` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close all views; Collapse all Databases.

---

### Browse-Node-Platform-01 — Platform: click each sub-node (P1, matrix)
**Goal:** each Platform sub-node clicks without errors.
**Type:** smoke + matrix
**Preconditions:** logged in as admin (`DATAGROK_LOGIN`); Context Panel expanded; Platform is expanded.
**Steps:** For each node from the list below:

**Groups (expand + click):**
- `Admin` (and all 6 inside: `Metrics`, `Reports`, `Service Logs`, `Test Manager`, `Test Track`, `Usage Analysis`)
- `Plugins`
- `Credentials` (and inside: `AWS`, `GCP`)
- `Functions` (and inside: `Queries`, `Scripts`, `OpenAPI`, `Api Samples`)
- `Predictive models` (and inside: `MLFlow`)
- `Dockers`
- `Sticky Meta` (and inside: `Types`, `Schemas`)

**Item leaves (single click):**
- `Users`, `Groups`, `Roles`, `Notebooks`, `MCP Servers`, `Repositories`, `Keys`, `Sync`, `Layouts`.

For each node:
1. Click the node in the tree. → If it is a group — it expands; if it is an item — the corresponding view opens.
2. `assertNoErrors(page)`. → No errors.
3. (For items) close the opened view before moving to the next one.

**Final check:** all ~20 sub-nodes clicked without errors.
**Selectors:** `treeGroupByName('Platform')`, `treeNodeByName(subnodeName)` (for each node in the list), `treeExpandArrow(groupName)`, `viewTabHandle(itemName)`, `VIEW_CLOSE`, `BALLOON_CONTAINER` — see `selectors.md` / `selectors.ts`.
**Postconditions/cleanup:** close all views; Collapse all Platform.

---

## Notes for automation (Playwright)

- **Base URL:** `https://dev.datagrok.ai/` → in the config's `baseURL`. Authorization — via `storageState`/global setup, not repeating the login in each test.
- **Case structure → code:** "Preconditions" = `beforeEach`/fixture; a step = action (`await ...`) + reaction check (`expect`); "Final check" = final `expect`; "Postconditions/cleanup" = `afterEach`/teardown.
- **Test tags:** use the "Type" field for tags — `@smoke`, `@functional`, `@regression`, `@negative`. The smoke set (all P1 smoke) — for running on every deploy.
- **Error detection:** a global `page.on('pageerror')` and `page.on('console')` listener (filter by `error`); plus checking for the absence of a datagrok error balloon (selector to be confirmed). Extract into a shared `assertNoErrors(page)` helper and call it at the end of every "no errors" case.
- **Creating data for preconditions:** create objects for Nav-06, Files-05, MyStuff-04, ModelHub-05 via API/second context in setup rather than by hand in the UI — this is more stable and faster.
- **Cleanup:** all cases with creation/renaming/share/saving clean up after themselves (see "Postconditions"); names with the `qa_autotest_<timestamp>` prefix + a shared teardown that deletes everything by prefix (a safeguard against failed tests).
- **Asynchrony:** the tree and lists load from the server — wait for the nodes/items to appear, do not `sleep`.
- **Timings:** for Apps-01, extract the 5 s threshold into a config constant so it is easy to change.
- **Selectors:** the real values are stored in `selectors.md`. On the first UI-recon run, capture stable `data-*`/roles for: the Browse/Favorites icons on the Sidebar; the tree nodes and expand/collapse arrows; the Top Menu (Home/Import file/Import text/Refresh/Collapse all/Locate); the context menu; the Favorites star and the Context Panel controls (Back/Forward/Collapse/Expand); the tabs and the pin button; the Hamburger (Split right/down); the mode switcher on the Status Bar; the global search bar; the filter panel and the search field. Under each case — a "Selectors" block with a reference to the logical names from `selectors.md`.
- **Accounts:** for permission cases (Tree-07, Platform-03, MyStuff-04) — a separate `storageState` under `DATAGROK_SHARING_LOGIN`. It is convenient to set up a second project in `playwright.config.ts` (`name: 'chromium-sharing'`) with its own `storageState`.
- **Matrix sets:** "open each application" and "open each demo file" — in the separate files `apps_matrix.test.ts` and `demo_files_matrix.test.ts` (tests parameterized over the result of the API request for the list of applications / files).

---

## Scope decisions

- **Spaces:** in Browse — only smoke (section 13). The rest — in a separate Spaces test set.
- **Model Hub / Compute:** smoke (section 14).
- **Global search / Context Panel / Split & modes:** added as sections 15–17 (basic coverage; deep scenarios — separately if needed).
- **Test data:** dev has demo files, demo DBs, and demo dashboards — the preconditions rely on them.

## Next steps

1. Walk through the P1 smoke cases by hand on dev and capture the real selectors → add a "Selectors" block to the cases.
2. Align the set with the team; enter it into Test Track if needed.
3. Start the port: first `@smoke` (P1), then `@regression` (bug cases).
