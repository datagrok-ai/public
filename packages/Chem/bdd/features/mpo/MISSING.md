# MPO qnum scoring — what the features cannot say yet

Written during the MPO qnum round of 16 Oct 2026 (public#4291, the manual checks of that day and the Slack
thread on MPO), updated on 20 Oct 2026. Covers `qnum-scoring.feature` and `score-dialog-cancel.feature`,
written with the existing steps only: no library, binding or code change. Items 7-8 were added on 28 Oct
2026 for `coverage-score.feature` (public#4403, checked on localhost with Chem and EDA from master); two more
questions of that round, answered on 29 Oct 2026, are under Won't do.

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
6. **API tests** (statistics `mapColumnDesirability` / `MpoCalculator`, Chem `mpoScoreByProfile`). Chem
   `mpo-tests.ts` (`MPO: qnum scoring`) has the three Score as options on a peaked curve, a bound wholly
   outside the curve, an inverted curve, a log10 curve with `0`, `<0`, `-2` and `>0`, the smallest exact
   value on a flat curve and the three missing-value strategies. Still missing:
   - the Score as options on rising and falling curves; a bound inside a flat part (the same in every option);
   - a log10 curve with `<-1` (empty score) and `>-1` (scored over the whole curve); `>0` with Limit value
     gives 0;
   - MPO Score by Profile on a qnum column uses the `qualifiers` the profile stores.
7. **Question: the column's name.** The option is "Coverage Score" (the name asked for in public#4403), the
   column it adds is `MPO <profile> coverage`. Is the difference intended?
8. **Package tests for the coverage count** (Chem `mpo-tests.ts` has one, `MPO: coverage`, over a plain, a
   categorical and a qnum column):
   - a value of 0 or below on a log10 curve is counted as measured although it gets no desirability;
   - a second `MpoCalculator.compute` with the coverage on, after one without it, counts every column (the
     cached desirability is not reused then);
   - `Chem:mpoCalculate` with `createCoverageColumn` returns the coverage column in its result table.

## Won't do

- **A check that the stored ADME Test is the shipped one**: the features run on a CI stand nobody edits by
  hand.
- **Every Score as option on every value, log10, MPO Score by Profile**: what the scoring computes, the API
  tests above.
- **The MPO Profiles app showing VALUES WITH < OR > only after the dataset is picked a second time**: not
  reproduced on localhost in three ways (create, then pick the dataset; edit a saved qnum profile; open it from
  the Browse tree, then pick the dataset). `qnum-scoring.feature` claims the first and the third.
- **The help page `mpo.md` saying nothing about qnum and Score as**: a help task.
- **An unchecked run removing a coverage column left by an earlier checked one**: OK with Coverage Score
  unchecked replaces the score column and keeps an `MPO <profile> coverage` column from an earlier run
  (`computeMpo` in `mpo/utils.ts` removes it only when the option is on), as intended (answer on public#4403,
  29 Oct 2026).
- **Coverage Score in the MPO context panel**: the option belongs to the MPO Score dialog only; the panel
  calls `computeMpo` with it off (`mpo/mpo-context-panel.ts`).
- **The coverage column as the marker size of a scatter plot** (the use asked for in public#4403): a scatter
  plot's Size property over an int column, nothing MPO adds.
