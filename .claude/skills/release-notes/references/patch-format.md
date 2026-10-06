# Patch release format (1.Y.Z, Z > 0, for example 1.27.12)

Patch releases use a **simplified format** — no intro paragraph, no Main updates,
no thematic sub-groups, no viewer sub-sections, no Packages section.

## Structure

```markdown
## <release-date> Datagrok <version> release

### Improvements:

* Component: description of improvement
* Another component: description
* Description without prefix (for platform-wide items)

### Fixed:

* Component now works correctly when ...
* Another positive outcome description
```

## Rules

1. **Two sections only**: `### Improvements:` and `### Fixed:` (with colons after heading).
2. **No intro paragraph** — go straight to `### Improvements:`.
3. **Flat bullet lists** — no nested sub-groups or viewer sub-sections.
4. **Component prefix** — use `Component:` prefix when the item is specific to a component
   (e.g., `Filter Panel:`, `JS API:`, `Scatterplot:`, `Grid:`). Omit prefix for
   platform-wide items (e.g., `Implemented EKS Pod Identity support`).
5. **Consolidate related commits** — multiple commits touching the same area should be
   merged into a single concise line (e.g., 3 scatterplot WebGPU commits →
   `Scatterplot: Improved WebGPU rendering performance`).
6. **Skip minor details** — patch notes should be concise. Omit items that are too
   granular or are minor behavioral adjustments within a larger feature.
7. **Omit sections if empty** — if there are no improvements, omit `### Improvements:`;
   if there are no fixes, omit `### Fixed:`.

## Reformulation rules (apply to patch and minor)

- Vary the leading verb (Introduced, Improved, Implemented, Exposed, etc.).
- Bug-fix phrasing: describe the positive outcome, not the problem
  ("now restores correctly", "no longer exposes").
- Bug descriptions under `### Fixed:` must NOT start with "Fixed".
