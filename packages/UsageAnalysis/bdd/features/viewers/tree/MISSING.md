# Tree edge cases — what the features cannot say yet

Written for `tree-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data shapes,
columns removed or renamed under the open viewer, layouts naming a missing column, property bounds,
closing the viewer. Nothing here asks for a library, binding or code change made by that feature; it
lists what the feature could not claim and why.

## Still open

### 1. Steps requested

None.

### 2. Questions for owners

1. **A column with more than 500 categories while the filter crosses 500 rows.** The category check
   runs only while at least 500 rows pass (`tree-viewer.ts:921-927`), so the viewer flips between
   its message and a 600-branch tree as the user filters, although the message says "fewer than 500
   unique categories" (`:931`). Should the check depend on the filter?
2. **Every row filtered out**: an empty tree with the root "All", or a message? (No error and the
   branches coming back are claimed.)
3. **The Size and Color columns removed.** The removal handler cleans the hierarchy only
   (`:629-633`); a removed Size or Color column then throws on every render (pinned, public#4423).
   Should Size and Color fall back to their defaults?
4. **A negative Label Rotate** (no minimum declared, `:130`), and values above the declared maximum
   typed in the panel: should the panel take them (see the trellis plot's MISSING.md)?

### 3. API or Dart tests to write

None.

## Won't do

Nothing.
