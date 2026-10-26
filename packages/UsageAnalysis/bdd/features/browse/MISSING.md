# Browse edge cases — what the features cannot say yet

Written for the edge-case features `address-edge-cases`, `file-preview-edge-cases` and `gallery-search-edge-cases`,
with no library, binding or code change. Each item below is a case those features leave out because the expected
behaviour has no source yet or the setup needs what the stand cannot give.

## Still open

### Questions for owners

1. **An empty (0-byte) file and an empty folder previewed.** No source says what the preview shows.
2. **An entity in the open tree deleted by another session, then clicked** (core#370: the tree shows a newly saved
   project only after Refresh). What should a click on a node whose entity is gone show — a "not found" message,
   or the node removed on the next Refresh only?
3. **A function address that names nothing** (`/func/<no such>`) says "Unable to find function" in a balloon and
   keeps the open view (claimed in `address-edge-cases`), and also logs the error with a stack trace to the
   console, while a project address that names nothing logs nothing. Is the logged error expected?

### Steps requested

4. **The address of an entity the user may read, inside a parent the user may not read (core#215).** The ticket has
   no description. The feature needs a second account that holds a grant on a child entity only (a project in a
   space it has no access to) and a step that opens an address as that account; the sharing steps sign in as the
   sharing user, but no step grants on the child alone without the dialog walk, and the expected outcome needs the
   ticket's own steps.

## Won't do

- **A file larger than the 10 MB preview limit** (`help/access/files/files.md`): the fixture would be written into
  the stand on every run; the limit is a reader's package test.
