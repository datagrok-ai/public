# PmxStudio changelog

## Unreleased

* The test (2026-11-10): each treated group is compared with the control with their own variances (Welch's t test
  on the per-animal slopes) and Holm's correction, in Plan (Crux `stats.tumor_study_design`, `test: welch`,
  `correction: holm`), Monitor (`stats.tumor_monitor`) and Analyze (Crux `stats.posthoc`; each difference's interval
  at 1 − α/k). With one variance for every group, an arm less spread than another lost 13 to 38 points of the power
  Plan gave it (three-group studies simulated with the engine: a small molecule beside an immunotherapy).
  * **Multiple comparisons**: `holm` (the default), `bonferroni` or `none`. Dunnett's test, which takes one
    variance for every group, is no longer a choice.
  * Analyze: beside each p-value, the one a common variance gives (`p_common_variance`: Dunnett's), and a note when
    the two disagree at α.
* Defaults (2026-11-10):
  * A new plan or study starts from a humane limit of 1500 mm³ and an endpoint of 1000 mm³ (2000 and 1500 before):
    the limit is under the 1.5 cm mean diameter that Workman et al. 2010 advise for therapy studies, a default study
    needs no more animals and ends 2 to 9 days sooner. A lab sets its own limit once, in the package's settings
    (`humaneLimit`); the plan says the limit's mean diameter and when it is above the guideline. Saved plans and
    studies keep their own limits.
  * Plan: the protocol says what the extra animals for losses give for every group at once, beside each group's
    own probability (`confidence_study` in the plan's table). The counts are unchanged.
  * Monitor: under the measurement sheet, how often a correct value is asked to be measured again and what that is
    for the study's animals in a session.
* Monitor, "measure again" (2026-11-10): one decision instead of two. The range to measure again outside, the
  expected range and the flag on an animal's last value all come from one forecast of its next measurement (Crux
  `stats.tumor_monitor`: `check_lower`/`check_upper`, `band_lower`/`band_upper`, `last_flag`), so the attention list
  and the sheet no longer disagree. Before, the flag ran at 3.7% and the range beside it at 2.6% for a nominal 1%.
  * The ranges' levels are the line's own when the engine priors carry them (`twin.conformal`), and the note under
    the sheet quotes how often a correct value fell outside on published controls left out of the fit; a line
    without that says so instead of quoting a number.
  * `pmxMonitorStudy`: `animals` has `last_flag` (-1 below, 0 inside, 1 above; empty on a first value) in place of
    `last_pit`, and `status` has `check_chance`.
* Review fixes (2026-11-10):
  * A study's details: saving them from the header no longer drops the notebook's notes, an ID given on saving, or a
    date of day 0 set from the calendar.
  * `.pmxstudy` files are saved by the columns' roles (a table whose columns are named otherwise saves and reopens),
    keep the table's other columns (caliper readings, cage, comments) and are checked when read; a damaged file of
    the same name no longer blocks saving.
  * Overview: a study opens from its file as it is now, the list is read again on coming back to it, and a second
    click does not open the study twice. In Plan, a plan saved under another name is another study, with its own ID.
  * Monitor: Save study and Analyze this study are in the ribbon; Propose groups asks before replacing groups and is
    in reach for a study with none yet; a planned review that falls on a day not measured is held on the next
    measurement; each arm's kill rate is kept in the plan and read from it (no design per forecast).
  * Analyze: the spider plot and the days to the endpoint follow every edit of the study table; a plan picked or
    taken off analyzes again; the evaluation day is a day with volumes; exports wait for the analysis in view.
  * Study table: pasted numbers are read or refused, never guessed ("1,234" is not taken as 1.234); the outlier
    rules leave excluded rows out, read a jump against what the animal's group did over the same days, and do not
    take a group's spread below a caliper's own error.
  * Charts: animal and group names from files are escaped in tooltips; the spider plot starts at the volume at
    randomization; a group hidden from a legend stays hidden.
  * `pmxMapPrior` failed for every tumor model; fixed.
  * Removed: the `MONITOR_READY` switch and its branches, unused exports, and the bdd features written for the
    interface before the redesign.
  * u2: a subscription or a control belongs to the control or the view that builds it; the package's CSS styles
    its own classes where it styled u2's; the package is recorded in u2's consumer baseline, and
    `grok check` passes.
* Monitor (`pmxMonitorStudy`, app route `/monitor`; needs Crux `stats.tumor_monitor` and `stats.tumor_allocate`): follow
  a running study, its groups hidden when the study is blinded. Today lists who needs attention and each animal's
  expected next value; Enrollment proposes balanced groups; Growth each group's forecast; Decisions flags a treated
  group with a low chance of a clear result, only on reviews fixed at the start, and logs what is decided. With an
  optional body weight column, Today lists animals near or over the weight-loss limit or not weighed, and Welfare
  shows Analyze's Tolerability on the weights so far. Followed studies are saved as `.pmxstudy` files; Plan has Follow
  this study.
* Overview: the app opens on its stages as cards (Plan, Monitor, Analyze) joined by the files that pass between them, the user's studies (saved plans and followed studies) with a demo portfolio, and Open a study: a plan, a followed study or a table of measurements, dropped anywhere or picked, opens in its stage. `/plan` opens Plan directly
* Plan: the animals a study needs as the assumed TGI goes (the Assumed TGI panel; through the function, the `tgiScenarios` input and the `scenarios` output)

## 0.0.1 (2026-10-26)

PMX Studio Preclinical, Plan and Analyze, on Crux (crate tgi) and the atlas of historical controls:

* Plan (`pmxPlanStudy`, app route `/plan`): de novo design of a xenograft or syngeneic efficacy study: end of study,
  kill rate per arm, analyzable animals per arm for the target power (assurance over the atlas's labs) with the
  multiplicity correction, animals to randomize and implant, power by group size, end of study and measurement
  schedule, the study calendar and the protocol. Plans saved as `.pmxplan` files (Save plan, Download; opened in Plan
  from Files).
* Analyze (`pmxAnalyzeStudy`, app route `/analyze`): the study analyzed as Plan pre-specifies it (per-animal slope of
  log volume, each treated group against the control), TGI and T/C, each treated animal against its untreated twin,
  the study against its plan, and tolerability from body weight (weight-loss grades, removals, NCI limits).
* The predicted and observed growth and the body weight are drawn by the package's own band chart (`band-chart.ts`,
  echarts): the area between a lower and an upper column, a line through it, one band per group.
* Run computes in the browser and saves nothing: the study's tables are not uploaded to the server.
* `pmxPriors`, `pmxMapPrior`: the atlas's control-growth priors and the MAP prior through Crux.
