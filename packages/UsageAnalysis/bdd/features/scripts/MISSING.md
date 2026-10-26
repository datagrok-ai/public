# Scripts edge cases — what the features cannot say yet

What `script-run-edge-cases.feature` leaves out, written for it: no library, binding or code change was made.
Each item names the case, why it is not a scenario yet, and what would make it one.

## Still open

### Questions for owners

1. **Closing the script editor with an unsaved edit.** `script_view.dart` has no close guard;
   `scripts-edit.feature` claims only that Save is disabled until the code changes. Should closing ask?
2. **The console call with an argument of the wrong type** (e.g. `BddScript("x")` for a dataframe input).
   What error should the console show?
3. **Renaming a script.** Scripts have no Rename... command (`script_meta.dart`); the name changes through the
   `name:` header in the editor. Is that the intended way, and should it refuse a name with no letter as the
   Rename dialog of queries and connections does? A project whose script is gone is claimed in
   `projects/projects-lifecycle-script.feature`.

### API or server tests to write

4. **Two scripts with one name.** `makeUniqueName` deliberately skips scripts, queries and projects
   (`repository_query.dart`), so a creation script resolved by name is ambiguous — the class of core#546.
   An ApiTests case once the owners decide.
5. **A header with extra whitespace or a malformed `input:` line** (core#200, core#604: the parsers drift apart).
   A parser unit test in ddt/xamgle, not a feature.

## Won't do

- **Interrupting a running script.** JavaScript scripts finish at once, and server-language scripts are not run
  in features (the bdd library's CLAUDE.md, "What never becomes a feature").
- **Empty, out-of-bounds and wrong-type scalar inputs of a script.** The run dialog of a script and of a query
  is the same input form; the cases are claimed once, in `queries/query-parameters-edge-cases.feature`.
