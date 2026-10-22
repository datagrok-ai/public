# Formula lines — what the features cannot say yet

Written 2026-10-22 when `formula-lines-dialog.feature` took scenarios 1 (GROK-17282) and 8
(GROK-19943) of TestTrack `PowerPack/formula-lines-dialog.md`. No library, binding or code change
was made.

## Still open

- **OK disabled in the Formula Lines dialog** (`formula-lines-dialog.md` scenario 1 step 3, GROK-17282).
  The dialog disables OK by adding a plain `disabled` class to the button
  (`src/dialogs/formula-lines.ts`, `onFormulaValidation`: `getButton('OK').classList.add('disabled')`),
  with no `aria-disabled` and no native `disabled`. The library's disabled reading
  (`libraries/bdd/src/runtime/assertions.ts`, `flagReadings`) knows `u2-input-disabled`, `d4-disabled`,
  `d4-filter-disabled`, `d4-menu-item-disabled` and `aria-disabled` only, so
  `OK button in "Formula Lines" dialog should be disabled` reads `false` while the button is greyed out.
  Wanted: the reading to recognise a dialog button's `disabled` class (or the dialog to set
  `aria-disabled`), then the scenario claims OK disabled after the bad formula and enabled after the fix.
