# PmxStudio changelog

## Unreleased

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
