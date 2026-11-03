# MPO qnum scoring — what the features cannot say yet

Written during the MPO qnum round of 16 Oct 2026 (public#4291, the manual checks of that day and the Slack
thread on MPO), updated on 20 Oct 2026. Covers `qnum-scoring.feature` and `score-dialog-cancel.feature`,
written with the existing steps only: no library, binding or code change. Coverage Score (public#4403) is claimed in
`calculate/mpo.feature`.

## Resolved

- **The smallest exact value of a qnum column scored 0 on a flat curve** (comment on public#4291). The
  default curve of a new profile starts at the column's minimum, read with the qualifier bits still in it,
  and the row's own value, read without them, fell just short of it. Fixed in `e51bdfc027`;
  `qnum-scoring.feature` claims row 1 scores like the rest, and SPGI's Competition assay scores 0.5 on all
  2486 non-empty rows.
- **Cancel in MPO Score kept the preview's score column** (public#4323). Fixed in `21599220c5`: the preview
  writes a hidden `~MPO <profile>` column and Cancel removes it; `score-dialog-cancel.feature` claims the
  table's own columns after Cancel.
- **Cancel in the curve settings kept the changed Score as** (public#4324). Fixed in `541f52dd10`;
  `qnum-scoring.feature` claims Score as is back and the score unchanged.
- **A bound wholly past the end of the curve kept the curve's end value** (asked on public#4291): on ADME
  Test's Lipophilicity, ending at 5 with 0.5, `>6` scored 0.5 while the exact `6` scored 0. Fixed in
  `360702be5d`: such a bound scores 0 like an exact value there, which `qnum-scoring.feature` claims.
- **Two features deleting each other's profiles.** The cleanup of `no MPO profile named {string} is on the
  server` deletes every profile the user saved during a scenario, so `profile-editor.feature` is `@serial`
  like `qnum-scoring.feature`. A narrower cleanup is still asked for below.

## Still open

1. **A property row by its name.** The `MPO property` kind (Chem `bindings/elements.ts`) matches the row's
   text, and in Design mode the name is an input, so a row is reached by its place. The order depends on the
   profile's source: the table's column order for a new profile; Caco2, Solubility, Lipophilicity for one read
   from the server; the saving order for one saved in the session. The features use Caco2, which is first in
   all three. Wanted: the kind matching the value of the property-name input, so that `"Settings" icon in
   "Caco2" MPO property` works.
2. **The stored Score as.** A step `the MPO profile {string} should store {string} for its {string} property`
   reading `properties[p].qualifiers` (`worst` / `limit` / `best`). The saved choice is proven through the UI
   only (reload, reopen, read the input); a server read would tell "not saved" from "not loaded".
3. **The log10 and Invert toggles.** Each is a `div.statistics-mpo-curve-btn` with no role and no
   `aria-label` (its tooltip is bound with `ui.tooltip.bind`), so no phrase reaches it. Wanted only if a
   "log10 + a bound" scenario is: an `aria-label` in the statistics library or two elements in Chem
   `elements.ts`.
4. **A cleanup of the feature's own profiles.** `deleteNewProfiles` (Chem `bindings/mpo.ts`) should delete only
   the profiles a feature names or whose save it observed. Until then every feature that saves MPO profiles
   has to be `@serial`, since the serial lane runs beside the parallel one (`playwright.config.ts`).
5. **OK's calculation.** A step that waits for the `Chem:mpoCalculate` call OK makes and fails if it threw.
   OK's errors are caught and shown as a balloon (`analysis/mpo.ts`, `runMpoCalculation`); the features claim
   the dialog closed and no error balloon, while the preview writes the same column.
6. **API and package tests**: listed in `libraries/bdd/API_TESTS_TO_WRITE.md` (the Chem MPO row): the Score as
   options on rising and falling curves and inside a flat part, log10 bounds, MPO Score by Profile on a qnum
   column, and the coverage count on a log10 curve, after a cached compute and in `Chem:mpoCalculate`.

## Won't do

- **A check that the stored ADME Test is the shipped one**: the features run on a CI stand nobody edits by
  hand.
- **Every Score as option on every value, log10, MPO Score by Profile**: what the scoring computes, the API
  tests above.
- **The MPO Profiles app showing VALUES WITH < OR > only after the dataset is picked a second time**: not
  reproduced on localhost in three ways (create, then pick the dataset; edit a saved qnum profile; open it from
  the Browse tree, then pick the dataset). `qnum-scoring.feature` claims the first and the third.
- **The help page `mpo.md` saying nothing about qnum and Score as**: a help task.
