# EDA Analyze — what the features cannot say yet

## Control comparisons dialog (closed-bugs round, 2026-11-05)

Both bugs were fixed in the core's tooltip and dropdown code; the dialog itself is in
`control-comparisons.feature`. Neither symptom has a reading a step could claim:

- **core#416: the tooltip of a red (invalid) Category field jumped left near the field's right edge.**
  With the dialog moved to the window's right edge and USUBJID picked as the category, the
  tooltip over the red field moved as the pointer neared the field's edge. A claim needs the
  tooltip's box at two pointer positions over the same field (`the tooltip should not have moved
  more than N px`), from the tooltip's position (`core/client/d4/lib/src/widgets/tooltip/`), and a
  step that moves a dialog to the window's edge. Neither exists.
- **core#418: the Control dropdown flashed black while expanding** with the red USUBJID category.
  A flash is a frame, not a state: no status reading or DOM state holds it after the list opens.
  A screenshot diff of the opening frames is the only check, and that is not a step the library
  offers.
