---
feature: spaces
target_layer: playwright
coverage_type: smoke
priority: p0
realizes: []
realized_as:
  - spaces-general-v2.test.ts
related_bugs: []
---

# Spaces — Playwright Test Reference

**File:** `playwright-tests/e2e/spaces/spaces.test.ts`
**Target:** `https://dev.datagrok.ai` (configured in `playwright-tests/.env`)
**Total tests:** 38 (38 active, 0 skipped)

All test actions are performed through the UI (clicks, keyboard, drag-and-drop, context menus, dialogs).
The JS API is used **only in `afterEach`** for cleanup, plus in two DnD tests for file creation
where no UI equivalent exists (documented per-test).

---

## Section 1 — Browse tree (3 tests)

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 1 | Spaces node is visible in the browse tree | The `.d4-tree-view-group-label` "Spaces" appears in the left panel within 20 s of loading |
| 2 | Clicking Spaces shows the spaces view with search bar | Single-clicking the Spaces tree node opens the spaces list view; `input[placeholder*="Search space"]` is visible |
| 3 | Right-clicking Spaces node shows "Create Space…" in context menu | Context menu on the Spaces node includes the "Create Space..." item |

---

## Section 2 — Create (5 tests)

Setup: all tests create spaces; `afterEach` deletes them via API.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 4 | Create root space via right-click on Spaces node | Right-click Spaces → "Create Space..." → dialog with "Create Space" title → fill name → OK → space card appears in list |
| 5 | Create child space via right-click on root space | Create root via UI; right-click root → "Create Child Space..." → fill child name → OK → child visible in content area |
| 6 | Create three nesting levels via UI | Root → child → grandchild, all created via UI; grandchild visible after navigating into child |
| 7 | Empty space name — OK button is disabled | Clearing the Name field: OK button loses `enabled` class; form cannot be submitted |
| 8 | Duplicate root space name — toast error shown | Creating a space with an already-existing name: toast "Root project with same name already exists"; only one entry in list |

---

## Section 3 — Read / Navigate (5 tests)

Setup: create space via `uiCreateRootSpace` helper; `afterEach` API delete.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 9 | Clicking a space shows its name and Details in the context panel | Single-click space card: entity name and "Details" / "Created" appear in right panel |
| 10 | Context panel shows all expected sections | Right panel has Details, Sharing, Activity, Content, Chats sections |
| 11 | Double-clicking a space opens it (breadcrumb updates) | Double-click opens the space view; breadcrumb/heading shows space name |
| 12 | Search filters visible spaces by name | Typing in the Search bar shows only matching spaces |
| 13 | Search with no match results in empty list | Searching for a nonsense string removes all space cards from the list |

---

## Section 4 — Rename (5 tests)

Setup: create space(s) via UI; `afterEach` API delete both ORIG and NEW names.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 14 | Rename root space — dialog pre-filled, list and tree update | Right-click → Rename... → "Rename project" dialog; input pre-filled with current name; new name appears, old name gone |
| 15 | Rename — browse tree node updates to new name | After rename, the `.d4-tree-view-group-label` shows the new name |
| 16 | Cancel rename — original name preserved | Clicking CANCEL in the rename dialog leaves the space name unchanged |
| 17 | Child space rename — dialog pre-filled with child name | Create child space; open Rename dialog on child; input pre-filled with child's name |
| 18 | Rename to existing space name — error shown, original name preserved | Rename a space to the name of another existing space → toast/error "same name / already exist"; original name still visible |

---

## Section 5 — Delete (3 tests)

Setup: create space via UI; `afterEach` API delete as safety net.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 19 | Delete root space — confirmation dialog, then space removed from list | Right-click → "Delete Space" → confirmation dialog "Are you sure?" → click DELETE → space card gone |
| 20 | Cancel delete — space is preserved | Open delete confirmation, click CANCEL → space card is still visible |
| 21 | Delete parent space — cascades to child (both gone) | Create parent + child; delete parent → both parent and child removed from list |

---

## Section 6 — Favorites (2 tests)

Setup: create space via UI; `afterEach` API delete.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 22 | Add to favorites — space appears under My stuff / Favorites | Right-click → "Add to favorites" → navigate to Favorites → space name visible |
| 23 | Remove from favorites — space no longer under Favorites | Add to favorites, then right-click → "Remove from favorites" → navigate to Favorites → space name absent |

---

## Section 7 — Share (1 test)

Setup: create space via UI; `afterEach` API delete.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 24 | Share dialog opens and shows user input + "Full access" | Right-click → "Share..." → dialog with user/group input, "Full access" label for current owner |

---

## Section 8 — Context menu completeness (1 test)

Setup: create space via UI; `afterEach` API delete.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 25 | Root space context menu has all expected items | Right-click on a root space: all five items present — Share..., Rename..., Delete Space, Create Child Space..., Add to favorites |

---

## Section 9 — Drag and Drop (11 active, 0 skipped)

**How DnD works in Datagrok:**
Dragging a file label from the Files view onto a space's tree node triggers a **"Move entity"** dialog
with a `select` dropdown (default: **Link**; options: Link / Copy / Move) and YES / CANCEL buttons.

**Setup:** all tests create `PW-DnD-Space` via UI; `afterEach` deletes it via API.
**File source:** `System.DemoFiles` (the built-in Demo connection at `/files/System.DemoFiles/`).
**Demo files used:** births.csv, cars.csv, cereal.csv, demog.csv, TSLA.csv, wells.csv, gapminder.csv.

**Important — Link vs Copy vs Move semantics:**
- **Link** adds the SAME entity to the space (shared ownership). Deleting the entity removes it from all spaces that reference it. Moving a linked entity from a space also removes it from all other spaces.
- **Copy** creates a NEW independent entity in the space. The original is unaffected regardless of what happens to the copy.
- **Move directly from Demo** removes the file from Demo permanently.
- Use **Copy** when setting up a disposable source entity for a Move test.

**Space-to-space drag timing note:** After clicking YES in the "Move entity" dialog for Link operations
between spaces, add `await page.waitForTimeout(600)` before clicking YES (let the dialog initialize)
and wait for the dialog to close (`waitFor({ state: 'hidden' })`) before navigating away.

**Important:** the `select` locator for the Link/Copy/Move dropdown is scoped using
`.filter({ has: page.locator('option[value="Link"]') })` to avoid matching other `<select>` elements on the page.

### Dialog structure

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 26 | Dialog has title, Link/Copy/Move options, Link is default (births.csv) | "Move entity" text visible; dropdown has value "Link"; all three option values exist |
| 27 | Dialog subtitle contains the target space name (cars.csv) | After drag, the space name appears somewhere in the dialog (subtitle "to [SPACE] project") |

### Link / Copy actions

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 28 | DnD Link — demog.csv appears in space | Drag + YES (Link default) → navigate to space → demog.csv visible |
| 29 | DnD Copy — TSLA.csv appears in space | Drag + select Copy + YES → TSLA.csv visible in space |
| 30 | DnD Move — births.csv moved from source space to target space | Create `PW-DnD-Source` space → **Copy** births.csv from Demo into it (Copy creates independent entity; Link+Move would remove it from Demo) → open source space → drag births.csv to `PW-DnD-Space` tree node with Move → births.csv visible in target space |

### Cancel

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 31 | DnD Cancel — wells.csv NOT added to space | Drag + CANCEL → navigate to space → wells.csv has count 0 |

### Multiple files

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 32 | DnD three files sequentially — births.csv, cars.csv, cereal.csv all in space | Three drags one after the other; all three visible in space content |

### Edge / atypical scenarios

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 33 | DnD same file twice — no error, at least one entry visible (gapminder.csv) | Drag the same file twice with Link; UI must not crash; at least one entry appears; no error toasts |
| 34 | DnD to space with existing content — new file added without replacing existing | Add demog.csv then beer.csv; both entries visible simultaneously |

### Delete-source effects

**Strategy:** Demo files cannot be deleted directly. Instead, **Copy** a Demo file into a source space
(creating an independent entity), then Link/Copy from source space to `PW-DnD-Space`, delete the source
entity (via right-click → Delete... → DELETE in source space), then verify target behavior.

**Context menu "Delete..."** opens a confirmation dialog ("This will delete the file for everyone and can't
be undone"). This permanently deletes the entity from the system, not just unlinks it from the space.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 35 | Delete source after Link — linked entry disappears from space | Copy wells.csv to source space → Link from source to `PW-DnD-Space` → verify in target → delete from source → verify entry **gone** from target (shared entity deleted) |
| 36 | Delete source after Copy — copied entry remains in space | Copy gapminder.csv to source space → Copy from source to `PW-DnD-Space` → delete from source → verify entry **still present** in target (copy is independent) |

---

## Section 10 — Forbidden (2 tests)

Setup: create space via UI; `afterEach` API delete.

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 37 | Duplicate root space name — toast error, only one entry in list | Create space; try again with same name → "same name" toast; exactly 1 entry in list |
| 38 | Empty name disables OK button in Create Space dialog | Open Create Space dialog; clear input → OK button loses `enabled` class; submit blocked |

---

## Helpers

| Helper | Purpose |
|--------|---------|
| `openSpacesView(page)` | Navigates to BASE URL, waits for Spaces tree node to be visible |
| `apiDeleteSpace(page, name)` | API cleanup: finds space by `friendlyName` via `grok.dapi.spaces` and deletes it |
| `rightClickSpacesTreeNode(page)` | Dispatches `contextmenu` on the Spaces tree node |
| `rightClickSpace(page, name)` | Dispatches `contextmenu` on a `.d4-tree-view-group-label:visible` node matching the name |
| `clickMenuItem(page, label)` | Clicks a `.d4-menu-item` by exact label text |
| `fillNameAndOk(page, name)` | Selects all text in the dialog input, fills new name, presses Enter |
| `uiCreateRootSpace(page, name)` | Full UI flow: right-click Spaces node → Create Space... → fill name → OK → wait for tree node |
| `dragFileToSpaceNode(page, fileName, spaceName)` | Navigates to `System.DemoFiles` via `page.goto`, scrolls file into view, drags to space tree node |
| `openSpaceContent(page, spaceName)` | Navigates to BASE, double-clicks space tree node, waits for `/s/` URL |

---

## Key CSS selectors

| Selector | Element |
|----------|---------|
| `.d4-tree-view-group-label` | Tree node labels (Spaces, Demo, PW-DnD-Space, …) |
| `.d4-tree-view-group-label:visible` | Visible-only tree nodes (used in `rightClickSpace` to skip hidden/clipped nodes) |
| `.d4-tree-view-node` | Tree node container (also has `d4-tree-drop` when drop target) |
| `.d4-menu-item` | Context menu items |
| `.d4-link-label label` | File/space cards in content area (shows full filename with extension) |
| `.ui-btn` | Platform buttons (OK, CANCEL, DELETE, YES) |
| `select` filtered by `option[value="Link"]` | Move / Link / Copy action dropdown (scoped to avoid other `<select>` elements) |
| `text=Move entity` | DnD dialog title |
| `input[placeholder*="User"]` | User/group input in Share panel |

---

## Coverage gaps (not automated)

| Scenario | Reason |
|----------|--------|
| DnD Move from Demo files directly | Move permanently removes the Demo file; workaround: first Copy to a source space, then Move from that space (used in test 30) |
| Drag root space into its own child | Requires drag; forbidden validation not verified |
| Drag from Dashboards to Space | Dashboards DnD not explored yet |
| Alt+drag (auto-link without dialog) | Modifier-key DnD not yet verified |
| Preview entities inside a space | Requires content; canvas-based preview not inspectable via DOM |
| Share space with another user | Needs second test account |
