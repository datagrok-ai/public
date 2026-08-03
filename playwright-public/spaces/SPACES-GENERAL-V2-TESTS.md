# spaces-general-v2.test.ts — Detailed Test Documentation

End-to-end Playwright tests for the **Spaces** feature in Datagrok. The suite covers the full lifecycle of spaces: creation, hierarchy, validation, rename, delete, favorites, sharing, search, drag-and-drop (copy/move), content view interactions, entity operations, and edge cases.

**Target URL:** `DATAGROK_URL` from `.env` (e.g. `https://dev.datagrok.ai`)  
**Browser:** Chromium (headless)  
**Workers:** 1 (sequential execution)  
**Default timeout:** 60s per test (some tests override to 120–180s)  
**Auth:** Pre-authenticated via `globalSetup` → stored in `e2e/.auth.json`  
**Sharing user:** `DATAGROK_SHARING_LOGIN` / `DATAGROK_SHARING_PASSWORD` from `.env`

---

## Helper Functions

### `openSpacesView(page)`
Navigates to the Datagrok home page and waits until the **Spaces** node is visible in the browse tree (up to 20s).

### `refreshSpacesTree(page)`
Clicks the sync icon (`i.grok-icon[name="icon-sync"]`) in the browse tree to refresh. Used after mutations to verify that refresh doesn't revert UI state.

### `apiDeleteSpace(page, name)`
Deletes a space by name using `grok.dapi.spaces` in the browser context. Used in `finally` blocks for cleanup.

### `rightClickSpacesTreeNode(page)`
Dispatches `contextmenu` on the root "Spaces" node. Uses `dispatchEvent` to bypass tooltip overlays.

### `rightClickSpace(page, name)`
Finds a space by name in the tree, dispatches `contextmenu` on its parent element.

### `clickMenuItem(page, label)`
Clicks a menu item matching text/regex, waits 600ms.

### `fillNameAndOk(page, name)`
In a dialog: clears input, types name character by character, presses Enter, waits 800ms.

### `uiCreateRootSpace(page, name)`
Right-clicks "Spaces" → "Create Space..." → fills name → asserts it appears in tree.

### `openSpaceContent(page, spaceName)`
Navigates home, double-clicks space in tree, waits for URL `/s/...` and 2s for content.

### `openSpaceViaTree(page, spaceName)`
Double-clicks space label in tree without full page reload. Used when already on Spaces page.

### `dragFileToSpaceNode(page, fileName, spaceName)`
1. Navigates to DemoFiles page
2. Checks file exists — if not, throws error with list of available files
3. Scrolls to file and drags onto space node in tree
4. Waits for "Move entity" dialog

### `addFileToSpaceViaLink(page, fileName, spaceName)` / `addFileToSpaceViaCopy(page, fileName, spaceName)`
Wraps `dragFileToSpaceNode` + clicks YES (Link) or switches to Copy + YES.

**Important:** All DemoFiles operations use **Copy** (not Link) to prevent deleting original files from DemoFiles.

### `rightClickItemInSpaceView(page, name)`
Dispatches `contextmenu` on an item in the space content view.

---

## Tests

### Test 1. Create root space and verify context menu

**Steps:**
1. Navigate to Datagrok, verify "Spaces" visible in tree
2. Click "Spaces" — stays visible
3. Right-click "Spaces" → click "Create Space..." directly → type `PW-Gen-Root-1` → Enter
4. Verify space appears in tree
5. Right-click space — verify menu has: Share, Rename, Delete, Create Child Space, favorites
6. **Cleanup:** delete via API

---

### Test 2. Validation: empty name, duplicate names

**Steps:**
- **Part A:** Empty name → OK button disabled
- **Part B:** Create `PW-Gen-Dup-2`, try duplicate → toast contains "name already exists", only one tree node
- **Part C:** Create parent + child, try duplicate child → toast contains "name already exists"

---

### Test 3. Hierarchy: child from tree, grandchild from view, navigate three levels

**Steps:**
1. Create root → create child **from tree** → verify child in tree (+ refresh, verify unchanged)
2. Double-click root → verify URL `/s/`, gallery grid visible, child card inside
3. Create grandchild **from content view** (right-click child card)
4. Navigate into child via dblclick → verify URL `/s/`, grandchild visible

---

### Test 4. Rename: pre-fill, cancel, success, duplicate error, child rename

**Steps:**
- **Part A:** Rename dialog pre-fills current name
- **Part B:** Cancel preserves original
- **Part C:** Successful rename — verify in tree + content view + right panel (+ refresh, verify unchanged)
- **Part D:** Duplicate name → error toast or original preserved
- **Part E:** Rename child space from content view — verify in view and tree (+ refresh)

---

### Test 5. Delete: cancel, selective, cascade

**Timeout:** 180s

**Steps:**
- **Part A:** Cancel delete → space preserved. Then delete from tree → verify gone from tree and view (+ refresh)
- **Part B:** Delete child1 from view → gone from view and tree, child2 survives (+ refresh)
- **Part C:** Cascade delete parent → parent and child2 gone from tree and view (+ refresh)

---

### Test 6. Favorites: add and remove

**Steps:**
1. Add to favorites → verify in browse tree Favorites section
2. Verify in sidebar favorites pane (star icon → `.grok-favorites-pane`)
3. Remove from favorites → verify gone from both locations

---

### Test 7. Share dialog: UI structure, share with permission, verify, delete removes access

**Timeout:** 180s

**Steps:**
- **Part A:** Dialog has User input, permission selector with "Full access" and "View and use"
- **Part B:** Select "View and use", type sharing user email, confirm OK
- **Part C:** API check — `grok.dapi.permissions.get()` confirms user in `view` permissions (not `edit`)
- **Part D:** Delete space → verify via API that space no longer exists
- **Part E:** Child space share dialog opens

---

### Test 8. Browse tree search

**Steps:**
1. Create space, click "Spaces", find search input
2. Type space name → card visible
3. Type non-matching string → card not visible

---

### Test 9. DnD dialog: options, default Link, Cancel, Link, Copy, duplicate add

**Timeout:** 180s

**Steps:**
- **Part A:** Drag `wells.csv` → dialog defaults to Link, has Copy/Move options. Cancel → file not in view or tree (+ refresh)
- **Part B:** Link `acidiq.csv` → visible in view. Click file → right panel shows info
- **Part C:** Copy `TSLA.csv` → both files visible. Click copied file → right panel shows info
- **Part D:** Duplicate add `acidiq.csv` → no error toast, exactly 1 instance

---

### Test 10. DnD Move: file moves to target, absent from source

**Timeout:** 180s

**Steps:**
1. Copy `beer.csv` to SRC, drag from SRC view to TGT tree with Move
2. `expect(dialogVisible).toBe(true)` — dialog must appear
3. Verify: gone from SRC view, present in TGT view
4. Click moved file → right panel shows info
5. Refresh → verify result unchanged in both spaces

---

### Test 11. Click file shows info in right panel

**Steps:**
1. Add `acidiq.csv` and `TSLA.csv` to space
2. Click `acidiq.csv` → panel shows "acidiq", has Details section
3. Click `TSLA.csv` → panel switches to "TSLA"

---

### Test 12. Click child space card shows info, panel switches between types

**Steps:**
1. Create two child spaces and add a file
2. Click child1 → panel shows child1 name, has Details
3. Click child2 → panel switches to child2
4. Click file → panel switches to file info
5. Click child1 again → panel switches back to space info

---

### Test 13. Double-click file opens table view; eye icon shows preview

**Timeout:** 120s

**Steps:**
- **Part A:** Double-click file → URL changes, ribbon and grid visible
- **Part B:** Go back to space, enable eye icon preview, click file → grid appears inside space view

---

### Test 14. Search inside space: exact, partial, no match, clear, child space

**Timeout:** 120s

**Steps:**
1. Add three files + create child space
2. **Exact match** "acidiq" → only acidiq visible
3. **Partial match** "aci" → acidiq visible
4. **No match** → all hidden
5. **Clear search** → all three files and child visible
6. **Child space search** → child visible, files hidden
7. **Verify DemoFiles intact** — `beer.csv` still exists in DemoFiles after test

---

### Test 15. Breadcrumb navigation: multi-level back, content preserved

**Timeout:** 120s

**Steps:**
1. Create child + grandchild (via view), add file to parent (via API)
2. Navigate: parent → child (dblclick card) → grandchild (dblclick card), verify URL `/s/` at each level
3. Navigate back: grandchild → child (via content view) → parent (via tree)
4. Verify: child card + file visible in parent after navigating back
5. Search still works after navigation

---

### Test 16. Delete file from space does not affect copy in another space

**Timeout:** 180s

**Steps:**
1. Copy `TSLA.csv` from DemoFiles to SRC, then DnD Copy from SRC to TGT (retry up to 3x)
2. Verify file in TGT
3. Delete from SRC → gone from SRC, still in TGT
4. Refresh → result unchanged

---

### Test 17. Entity ops: context menus, rename, cancel rename, delete, cancel delete

**Steps:**
1. Create space with `TSLA.csv` (Copy), `acidiq.csv` (Copy), child space
2. **File context menu:** Open, Rename, Delete present
3. **Rename file:** `TSLA.csv` → `PW-Gen-tsla-renamed` (+ refresh)
4. **Cancel rename:** try rename `acidiq.csv`, cancel → unchanged
5. **Cancel delete:** try delete `acidiq.csv`, cancel → unchanged
6. **Delete file:** delete `PW-Gen-tsla-renamed` → gone, acidiq still there (+ refresh)
7. **Child context menu:** Rename, Delete, Share present
8. **Rename child:** from view, verify in view and tree (+ refresh)
9. **No "Duplicate"** in parent space context menu

---

### Test 18. Edge cases: name with spaces, circular drag, duplicate child name

**Steps:**
- **Part A:** Create space `PW Gen Name With Spaces` → visible in tree
- **Part B:** Create parent + child, try duplicate child → error toast
- **Part C:** Drag parent onto child → app doesn't crash, parent still visible

---

## Data Isolation

- All space names use `PW-Gen-*` prefix with test number suffix
- All tests clean up via `try/finally` + `apiDeleteSpace`
- Tests 9, 10, 16 also pre-delete spaces at start
- All DemoFiles operations use **Copy** (not Link) to protect originals

## Key Design Decisions

- **`refreshSpacesTree` pattern:** Check result WITHOUT refresh first, then refresh and verify unchanged — catches bugs where refresh reverts state
- **DnD retry loops:** Up to 3 attempts for DnD operations that can fail due to tooltip overlays
- **DemoFiles file check:** `dragFileToSpaceNode` throws descriptive error if file not found, listing available files
- **No Link from DemoFiles:** All `addFileToSpaceViaLink` replaced with `addFileToSpaceViaCopy` to prevent deleting DemoFiles originals

## Selector Patterns

| Selector | What it matches |
|----------|----------------|
| `.d4-tree-view-group-label` | Labels in the browse tree |
| `.d4-menu-item` | Context menu items |
| `.d4-dialog` | Modal dialogs |
| `.d4-link-label` / `.d4-link-label label` | Item cards in content views |
| `.d4-toast`, `.d4-balloon` | Toast/balloon notifications |
| `.grok-prop-panel`, `.d4-info-bar` | Right panel (context/property panel) |
| `.d4-accordion-pane-header` | Accordion section headers (Details, etc.) |
| `.ui-btn` | Dialog buttons (OK, CANCEL, DELETE, YES) |
| `#elementContent .grok-gallery-grid` | Gallery grid in content view |
| `.grok-favorites-pane` | Sidebar favorites panel |
| `i.grok-icon[name="icon-sync"]` | Refresh/sync icon |
| `i.grok-icon[name="icon-eye"]` | Preview eye icon |
| `input[placeholder*="Search"]` | Search input |
| `input[placeholder*="User"]` | User input in share dialog |
