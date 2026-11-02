# Markup edge cases — what the features cannot say yet

Written for the markup edge cases, now an `@edge-cases` scenario of its own, after the master scenarios of `markup-person-paths.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Questions for owners**

None. Html mode content from a shared layout or project runs its scripts (`NodeTreeSanitizer.trusted`): core#1239, no scenario (a script run is not observable from a step).

**API or Dart tests to write**

1. Dart regTest: two `{{…}}` expressions on one line (`{{a}} and {{b}}`) are captured as one by the greedy `{{(.*)}}`.
2. ApiTests: a cell value holding HTML (`<b>bold</b> <i>it</i>`) substituted into `${NAME}` is shown as written text, with no bold element rendered (`d4/CLAUDE.md` rule 8, `htmlEscape` in `markup_viewer_core.dart`).
3. ApiTests: a column removed while the content references it leaves the reference as written (`Age: 26 Sex: ${SEX}`), no error.
4. ApiTests: once every row of the table is deleted, `#{t.rowCount}` reads 0, the viewer has no current row and reports no error.
5. ApiTests: a referenced column renamed (AGE to YEARS) is followed by the content (`Age: ${YEARS}`) and its values are still substituted: a known failure, core#1232, still open.

Moved out of BDD (no user gesture); requested from ApiTests.

## Won't do

- Long content and numeric bounds: the viewer has no numeric property and long text is plain DOM text.
