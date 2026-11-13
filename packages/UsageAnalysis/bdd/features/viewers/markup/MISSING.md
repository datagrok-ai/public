# Markup — what the features cannot say yet

Written in the BDD PR round of 8 Nov 2026.

## Resolved

- Two `{{…}}` expressions on one line were captured as one by the greedy `{{(.*)}}` and rendered as nothing: Dart
  `Markup substitutes two {{...}} expressions on one line one by one`.

## Still open

- **core#1239** — Html mode renders content from a layout or a project as trusted HTML (`NodeTreeSanitizer.trusted`),
  so its scripts run: the owner's call; no scenario, as a script run is not observable from a step.

## Won't do

- Edit content… dismissed with CANCEL or Escape: every dialog writes on OK only.
- A cell value holding HTML: Dart `Data shows as text, not markup…`.
- A column removed or renamed, every row deleted: `markup-row-binding.feature` "A reference to a column that does not
  exist is left exactly as written", Dart `Markup viewer content follows a column rename`, `Data Manipulation | Markup`.
