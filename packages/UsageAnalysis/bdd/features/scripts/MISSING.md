# Scripts edge cases — what the features cannot say yet

What `script-run-edge-cases.feature` leaves out, written for it: no library, binding or code change was made.
Each item names the case, why it is not a scenario yet, and what would make it one.

## Still open

### Waiting for core#456

1. **Closing the script editor with an unsaved edit**: the confirmation is requested in core#456.

### Steps requested

2. **A script renamed through its `name:` header** (decided 2026-10-27: that is the way to rename a script, as
   `script_meta.dart` has no Rename... command, and a name with no letter is refused as in the Rename dialog of
   queries and connections). Saving a script does not check the name today (`isValidEntityName` is used only by
   `DbEntityMeta.getRenameDialog` and the project Save dialog). No step edits one line of the editor: a step
   `user replaces {string} with {string} in {element}` would let a scenario turn `name: <script>` into
   `name: 2026` and claim the refusal.

### API or server tests to write

3. **Two scripts with one name.** `makeUniqueName` deliberately skips scripts, queries and projects
   (`repository_query.dart`), so a creation script resolved by name is ambiguous — the class of core#546.
   An ApiTests case once the owners decide.
4. **A header with extra whitespace or a malformed `input:` line** (core#200, core#604: the parsers drift apart).
   A parser unit test in ddt/xamgle, not a feature.

## Won't do

- **Interrupting a running script.** JavaScript scripts finish at once, and server-language scripts are not run
  in features (the bdd library's CLAUDE.md, "What never becomes a feature").
- **Empty, out-of-bounds and wrong-type scalar inputs of a script.** The run dialog of a script and of a query
  is the same input form; the cases are claimed once, in `queries/query-parameters-edge-cases.feature`.
