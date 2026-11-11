# PmxStudio changelog

## Unreleased

* Plan, second pass (2026-11-11): the animals per group and the end of study leave the form: a click on the
  Analyzable or End of study card sets them, the card says "set", and another click changes them or goes back to the
  computed value. Multiple comparisons moves to Advanced (Holm by default, none needed with one arm; the protocol names
  it). The calendar's instructions show while the pointer is on it, the first three times. When the controls may level
  off below the endpoint and no lower endpoint avoids it (syngeneic lines), the warning says to plan to end on the
  study's day instead of suggesting a lower endpoint.
* Plan's arms (2026-11-11): each arm's name is its card's title, and under it, in grey, what Plan reads in it (drugs,
  dose, schedule, route) and the kind of treatment it takes, a click to change it. A name Plan cannot place says so;
  in immunocompetent mice (or before a line is chosen) the card asks whether the compound works on the tumor directly
  or through the immune system, and assumes the latter, which needs more animals, until answered. A combination with
  a drug Plan cannot place asks the same for that drug. TGI or T/C is set once for the study (`effectAs`, Advanced;
  TGI for a new plan). When the controls may level off below the endpoint, the warning offers a lower endpoint that
  avoids it and applies it on a click.
* Plan (2026-11-11): a line with several equally good fits (Panc02, PDX-PDAC) is designed in a fraction of the time
  (two arms: 2.7 s to 0.6 s on Panc02), the other fits sharing one design's simulated studies; the group size is
  stated as estimated by simulation, accurate to about one animal (protocol and the Analyzable header's tooltip).
* Atlas (2026-11-11): the engine priors are the atlas sealed that day (`files/engine_priors_v3.json`, sha256
  `e2c069fe…d605`), in place of the export of 2026-11-06: the curation notes applied and three exclusions of
  non-volume or misaligned rows (Panc02, Hepa1-6), the measurement error recalibrated (4T1's tails), the digital
  twin refitted on the final data, each line's `twin.conformal` and `twin.next_ranges` (Monitor), and
  `calibration.identification` on the lines whose spreads the controls do not pin down. `files/SOURCES.md` lists
  every paper behind the twelve lines (219: each with a control in a line's pool or a treated arm behind an
  archetype); two papers were attributed to the wrong PDX line and one had no DOI. Known issue: a line's
  non-progressing share (`take`) is 6 to 12% low, relatively (0.03 to 0.1 points), from Crux's binomial
  `stats.map_prior` on tables with few events; no group size depends on it.
* The test (2026-11-10): each treated group is compared with the control with their own variances (Welch's t test
  on the per-animal slopes) and Holm's correction, in Plan (Crux `stats.tumor_study_design`, `test: welch`,
  `correction: holm`), Monitor (`stats.tumor_monitor`) and Analyze (Crux `stats.posthoc`; each difference's interval
  at 1 − α/k). With one variance for every group, an arm less spread than another lost 13 to 38 points of the power
  Plan gave it (three-group studies simulated with the engine: a small molecule beside an immunotherapy).
  * **Multiple comparisons**: `holm` (the default), `bonferroni` or `none`. Dunnett's test, which takes one
    variance for every group, is no longer a choice.
  * Analyze: beside each p-value, the one a common variance gives (`p_common_variance`: Dunnett's), and a note when
    the two disagree at α.
* Plan, lines the atlas does not pin down (2026-11-11): where the published controls are described equally well by
  several fits of how much animals and studies differ (the engine priors' `calibration.identification`: today
  Panc02 and PDX-PDAC), an arm's power is the average over those fits, as it already is over the atlas's labs, and
  the group is sized for it; a note gives what each fit needs alone when they differ, and says nothing when they
  agree. The power by measurements per week carries what the other fits move the group by at the study's own
  schedule (they are not designed again per schedule), and the power by last day and the smallest detectable
  effect stay the prior's own fit's. `pmxPlanStudy`: `armResults` has `n_fits`, what each fit needs alone (the
  prior's own first), and **Equally good fits** (Advanced) sizes for the atlas's best fit alone when asked.
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
  * The note under the sheet quotes how often a correct value fell outside the range shown before measuring, on the
    atlas's published controls (the engine priors' `twin.next_ranges`, an average per study), and on which side.
    The ranges' levels are Crux's own: the levels fitted per line (`twin.conformal`) are those of the flag on a last
    value, whose forecast is built after the session, and Crux takes one pair of levels for both; applied to the
    ranges shown before measuring they would leave a correct value outside twice as often as asked.
  * `pmxMonitorStudy`: `animals` has `last_flag` (-1 below, 0 inside, 1 above; empty on a first value) in place of
    `last_pit`, and `status` has `check_chance` and `check_above` (the share of those values above the range).
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
