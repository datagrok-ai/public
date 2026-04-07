# Spaces — Missing Test Coverage

Analysis based on `spaces-general-v2.test.ts` (18 tests, all passing) and [wiki](https://datagrok.ai/help/datagrok/concepts/project/space).

---

## Drag and Drop

| Scenario | Status |
|----------|--------|
| DnD from local storage (single file) → space | ❌ Missing |
| DnD from local storage (group of files) → space | ❌ Missing |
| DnD space→space via **Browse Tree** (Link / Copy / Move) | ❌ Missing |
| DnD **group of files** from space view → another space | ❌ Missing |
| DnD from **Dashboards**: project with linked tables → space | ❌ Missing |
| DnD from **Files**: folder → space | ❌ Missing |
| DnD from **Packages** → space | ❌ Missing |

---

## Read / Navigate

| Scenario | Status |
|----------|--------|
| Preview **project** in space | ❌ Missing |
| **Open project** from space | ❌ Missing |
| Preview linked file / linked project / linked function | ❌ Missing |

---

## Share

| Scenario | Status |
|----------|--------|
| Share **linked entity** | ❌ Missing |
| Cross-user: shared space visible in other user's Browse tree | ❌ Missing (verified manually, not automated) |
| **Child space inherits permissions** from parent | ❌ Missing |

---

## Delete

| Scenario | Status |
|----------|--------|
| Delete **linked project** / **linked function** from space | ❌ Missing |
| Delete **linked space** (space linked inside another space) | ❌ Missing |

---

## Forbidden / Edge Cases

| Scenario | Status |
|----------|--------|
| Drag child space to root Spaces level | ❌ Missing |
| Drag a function from Packages → space (should be blocked) | ❌ Missing |

---

## Rename

| Scenario | Status |
|----------|--------|
| Rename entity from Browse Tree | ❌ Missing — tree file nodes only render after space content is loaded and may appear hidden; unreliable in automation |
| Delete entity from Browse Tree | ❌ Missing — same limitation as Rename |


---

## Notes on Demo File Safety

All tests use `addFileToSpaceViaCopy` (not Link) when working with DemoFiles to prevent permanently removing original files. `dragFileToSpaceNode` helper throws a descriptive error with available file list if a DemoFiles file is missing.

---

## Priority Summary

| Priority | Area |
|----------|------|
| 🔴 High | Link/Clone semantics (read-only check, edit propagation) |
| 🟡 Medium | Open/preview project in space, preview linked entities |
| 🟡 Medium | Forbidden — drag child→root, packages→space |
| 🟡 Medium | Cross-user sharing automation, child space permission inheritance |
| 🟢 Low | Advanced metadata search, multi-select DnD, local file upload DnD |
