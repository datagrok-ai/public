# Browse edge cases — what the features cannot say yet

Written for the edge-case features `address-edge-cases`, `file-preview-edge-cases` and `gallery-search-edge-cases`,
with no library, binding or code change. Each item below is a case those features leave out because the expected
behaviour has no source yet or the setup needs what the stand cannot give.

## Still open

### Decided and written

An empty (0-byte) file is previewed with no error (decided 2026-10-27): `file-preview-edge-cases`, a .csv and a
.txt.

### Steps requested

1. **An empty folder previewed** (decided: no error). No step makes a folder without a file in it; a step
   `a folder {string} is on the server` (files API, deleted at feature end) would let `file-preview-edge-cases`
   click it.
2. **The address of an entity the user may read, inside a parent the user may not read (core#215).** The ticket has
   no description. The feature needs a second account that holds a grant on a child entity only (a project in a
   space it has no access to) and a step that opens an address as that account; the sharing steps sign in as the
   sharing user, but no step grants on the child alone without the dialog walk, and the expected outcome needs the
   ticket's own steps.

## Won't do

- **A file larger than the 10 MB preview limit** (`help/access/files/files.md`): the fixture would be written into
  the stand on every run; the limit is a reader's package test.
- **An entity deleted by another session, then clicked in the open tree**: it needs two sessions at once, which CI
  cannot hold.
