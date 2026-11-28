# PmxStudio changelog

## 0.1.3 (unreleased)

* Analyze: one test, one interval. Each treated group's result and its uncertainty come only from the growth rates:
  **Growth rate vs control: grows 4.1%/day slower (95% CI 2.0–6.3%), p = 0.003** (the interval adjusted as the
  p-value, its level said). TGI and T/C are the evaluation day's numbers, without an interval, in the groups table,
  the group's panel, Plan vs observed, the exposure chart, the results in words, the CTD table and the exports; the
  note that explained a T/C interval at odds with the result is gone.
* Monitor: **Enter today's values** is always there. On a day off the plan's calendar it asks first; with no plan it
  opens at once; with no dated day 0 it asks for the day number. A day entered off the calendar is a measured day like
  any other.

* One view: Plan, Monitor, Analyze and Compare replace each other in one tab, and coming back to a stage finds it as
  it was. A stage with unsaved changes stays open. Back and Forward walk the pages inside it, one press a step; a
  link to a page, or the app opened again, moves the open view instead of opening a second one.
* The view is always named **PMX Studio** (its title and the browser tab); the study is named in the page's header.
  The ribbon shows a home icon, back to the Overview, then the steps Plan — Monitor — Analyze. A stage with nothing
  to show yet is greyed, and its tooltip says why.
* Compare is no longer a step: it opens from the Overview (tick two studies or more, then Compare; the ticks are kept)
  and from Analyze (Study > Compare with other studies…). Its ribbon reads ⌂ Compare · N studies.
* The view's menu is Datagrok's own: **File** (Open… Ctrl+O, Save Ctrl+S, Save as…, Export, Documents), View, Study
  and Help, in place of the Open and Save buttons. An item that is not available yet says why on hover.
* Follow this study and Analyze this study are outlined, so they no longer look like the page's main action.
* Each panel's explanation is behind an (i). Notes that repeated a button or a field are gone.
* Panels fold to one line that says what they hold. Plan folds the calendar, Assumed TGI, What if and the protocol;
  Monitor the calendar and the data checks; Analyze time to endpoint, exposure, growth rates, notes and the spider;
  Compare the exposure and the atlas controls. Your folds are remembered, and Reset restores the defaults.
* **Customize this page** works in place, from the toolbox or View: drag a panel by its header (a long press
  on touch) and the others slide out of the way; release it anywhere and it settles into the slot shown. Beside a narrower panel with room left, it
  takes that room while it hovers there and keeps it when released; moved away, it springs back to its own width. Drag its right
  edge and the width snaps to a quarter, a third, a half, two thirds or the row. × hides a panel. Escape cancels a
  drag, Ctrl+Z undoes the last change. With reduced motion, everything moves at once. A hidden panel stays listed in
  the toolbox with + show, and in a line at the end of the page.
* The toolbox is the page's index: a click opens the panel's tab, unfolds the panel and scrolls to it.
* Linear · Log is plain text in each growth chart's header, the one in use underlined.
* Monitor and Analyze have two cards each; the rest is in the status line. Monitor opens on Needs attention.
* Plan: To randomize is under the Implant card. A click on a group opens its doses, power and Match a total… in the
  context panel.
* Analyze: a click on an animal shows its curve against its twin in the context panel.
* Compare: one summary line in place of the cards. The treated groups are listed under their study, each in its
  colour, the numbers aligned and sortable, where a number comes from on hover. Dose and effect: a colour per drug, a
  mark per study and the key above the plot (no labels over the points); it says what to add when it has nothing to
  draw.
* More space between panels.
* PMX Studio is listed under Preclinical in Browse.
* Saved plans and studies (schema 1.5.0): the document is stored as JSON; the owner is the Datagrok user, and whose
  a copy is, a user too (or, for someone else's file brought in, the name the file gives); day 0 is a date and the
  study days whole numbers. Study details' Owner is now **Run by**, free text, empty unless you fill it in. Project is
  gone: Datagrok has no way to put a saved row in a Space, and an earlier project is kept as a tag. Rows saved with
  0.1.2 are converted when the package is deployed (`databases/pmx/0001_schema_1_5_0.sql`).

## 0.1.2 (2026-11-26)

* A plan whose controls never reach the endpoint within the study (Panc02 at 2000 mm³) saves and opens again; when a
  plan or study of the app cannot be saved, the message says why.

* Analyze: NCA's exposure has a panel of its own, **Exposure and effect**, out of the groups panel; the groups table
  gains a **× patients' exposure** column when the patients' exposure is given.
* Numbers in exposure, AUC, Cmax, concentration, dose and multiple tooltips: 2–3 significant figures, no float noise
  ("0.30000001192092896" is now 0.3), thousands separators.
* Animals by day: one line under the growth chart; each point's tooltip says "n of N animals", and a stretch with fewer
  than 3 animals is dashed. The full table is in the group's context panel, the report and **Open as table**.
* A toolbox in Plan, Monitor, Analyze and Compare: the view menu's items by section (Study, Export, Documents, View,
  Help). Analyze exports the study file (`.pmxstudy`) and the calendar (`.ics`); Monitor's Propose groups… and the
  randomization sheet are in its Study section, Enrollment before the measurement sheet.
* Right-click the platform's way: animals, groups and figures are Datagrok objects with registered actions, the same
  in the right-click menu and in the context panel's Actions; Include in the analysis is offered for an excluded
  animal. The custom study menus are removed.
* Compare v2: panels you can lay out, with summary cards, every group's growth (Linear, Log, Relative V/V0), dose and
  effect across studies, TGI against the patients' exposure, each control against the atlas for its model (a table and
  one sentence when a control drifts), and the treated groups with their result.
* One Crux module, `pmx` (crux-core main 2623118b4), in place of `nonclinical` and `stats`: 5.6 MB of wasm instead
  of 21 MB, the same results.

## 0.1.1 (2026-11-26)

* Browse stays open: opening PMX Studio, or coming back to its Overview, no longer closes Datagrok's Browse tree,
  toolbox or help. Pin the PMX Studio view to give it the whole screen: Browse closes while the pinned view is current
  and comes back when you leave it or unpin it.
* One dialog at a time: clicking a button again while its dialog is open (Open a study, share, delete, a model's
  history, the calendar, ARRIVE, patients' exposure, the layout, the study's details and the others) brings the open
  dialog to the front instead of opening a second one.
* Open… in Plan, Monitor and Analyze: the same Open a study dialog as on the Overview, for that stage; files dropped on
  a stage's page open there too. In Plan, a table can be the growth history of a model outside the atlas. A file of
  another stage opens in its own stage, and the dialog says so.
* A clearer summary once a file is read: the study, its groups, animals and days, the files it came from, what was read
  and how, warnings for anything not read, and one main button with the other stages beside it.
* Demo data one click away: **Try with demo data** at the foot of each stage's card on the Overview, and, while Your
  studies is empty, **Or compare two demo studies** (CT26 and MC38 trametinib, saved as yours and opened in Compare).
* The Overview's toolbox: Open a study, New plan, Compare saved studies, the demo data and the help.

## Unreleased

* The priors file holds aggregates per tumor model only (2026-11-16): the per-paper table (`studies`) is no longer
  distributed, and `pmxMapPrior`, which recomputed the MAP prior from it, is removed; the MAP prior shipped in the file
  is the one used. The example study downloads as `pmx-example-mc38-yu2024-ccby.csv`, its source in its name.

* Groups balanced on body weight too (2026-11-13): when the study table has a body weight column, Propose groups
  balances on tumor volume first, then on body weight, and shows each group's mean and spread; if some animals lack a
  weight it says how many and balances on tumor volume only, and with complete weights it says body weight is not
  used yet while the engine lacks it.

* One reader for tables of measurements (2026-11-13): Open a study reads one row per measurement, a column per day or
  a column per animal, and says what it read before opening: the groups, the days, the units (converted, or assumed
  and said: mm, mm³, g), day 0 (a randomization or treatment-start date column, else the first measurement) and the
  volumes (as given; computed from length and width, L×W²/2 unless another formula is chosen; or the lab's formula
  found from its own volumes). Dates and week numbers become days. Every number is read or listed as not read, with
  the reason; what it cannot tell (two volume columns, day or month first, an implausible unit) it asks. A blinded
  study's summary names no group. Where a study's volumes come from is saved with it.
  Excel (.xlsx) files open too, through the PowerPack package: with several sheets, the one that reads as a study, or
  the user's pick.

* Combinations in Plan (2026-11-13): a group named like the sum of two others ("A", "B", "A + B") is taken as their
  combination, without asking (a tick box on the arm turns it off). Plan adds a note with how many animals per group
  would show the combination beats each drug alone (the group sizes stay as planned) and a clause on the synergy it
  could detect. In Analyze, Combination synergy… is an exploratory test (the plan did not set it): synergy at the plan's
  α, in plain words.

* Analyze and Monitor check the study's end against the same widened range Plan prints for MC38, CT26 and B16F10, so a
  study ending inside Plan's range is never called earlier or later than predicted.

* Plain messages from Open a study to Monitor (2026-11-12): what is missing, empty or failed says what happened and
  what to do; no technical text reaches the user (it goes to the console), and no failure is dropped silently (a
  saved-plans list that could not be read no longer says there are none). The atlas file is downloaded once.

* Provenance (2026-11-12): a saved plan or study records the atlas it was made with (a fingerprint of
  files/engine_priors_v3.json) and the PMX Studio version, an optional field: older files open as before. Opening a
  plan saved with another atlas says so in one line (Plan recomputes it with the current one; Monitor keeps the plan's effects and end as saved,
  its forecasts use the current one). The PDF reports name the app version and the atlas, and the plan's when it differs.

* Ethics application text (2026-11-12): beside the study protocol's Copy, a text for an animal-ethics application:
  why this number of animals, why not fewer, the losses expected and how they are covered, and reduction, with the
  plan's own numbers only and a bracketed line for the researcher to justify the expected effect.

* Plan's bands as validated (2026-11-12): in MC38, CT26 and B16F10 the predicted tumor growth range and the range
  of the end of study ("in 8 of 10 studies") are drawn as wide as the validation on published studies left out of
  the atlas found them (files/band_calibration_v1.json, with its provenance): the growth range is where about 9 in
  10 measurements of an animal are expected. The animals per group and the power are unchanged. Other lines keep
  the engine's ranges; a plan saved before keeps the ranges it was saved with.

* A study keeps its tumor model (2026-11-12): saved without a plan, it opens again with the model chosen in Monitor
  and the Overview names it (an optional field: older files read as before; with a plan, the plan's model is read).
  A blinded study names no arm or compound outside itself until its groups are revealed: not in the Overview, its
  search or the file preview, even when the blinding is in its plan's details only.

* Fixes (2026-11-12): in Plan, the measurement schedule comparison no longer spins without end when it is reached by
  scrolling; it fills in about a second. In Analyze's Each animal, the legends of the waterfall and the spider plot stay on one line,
  paged, instead of covering the bars and lines when the groups' names are long. In the Overview, a study opened
  from a table, with no plan or details, is listed by its groups (not while they are blinded).

* Exposure (2026-11-12): with the NCA package installed, Analyze's Advanced fields take a table of plasma
  concentrations (picked, or opened while the field is empty) and show each group's AUClast and Cmax from NCA, mean
  and CV%, beside its effect. A group is matched by its own name in the table, else (one drug) by drug, dose and route
  read in its name; nothing is guessed or scaled across doses or schedules, and a group NCA could not compute or
  match shows "—" and says why. With 3 doses or more, a chart of exposure against TGI, one point per group, no line
  and no test. Nothing is saved in the study, and the analysis is unchanged.

* Files (2026-11-12): a plan or a study selected in Datagrok's Files shows what it holds, read from the file alone:
  a plan its model, end, groups and animals; a study its animals, groups (not named while blinded), days measured,
  marked values and plan. A button opens it in Plan or Monitor; a damaged file or one from a newer version says so.

* An example study (2026-11-12): the empty Monitor and Analyze offer a published one to download, MC38 with LD4172,
  anti-PD-1 and the two combined (Yu et al., Nat Commun 2024, CC BY 4.0; source and changes in
  `files/example/README.md`).

* The Overview's list (2026-11-12): a folder that fails to list keeps the list shown and says so, instead of an empty
  list ("No studies yet"); a file in the plans or studies folder that is no PMX Studio file is read once, not on
  every Overview (one of ours from a newer version is read each time: an update reads it).

* The app starts clean (2026-11-11): the Overview lists only the user's own studies ("No studies yet" until there is
  one) and computes nothing when it opens; the demo studies are gone (the package's tests keep the simulated one, in
  tests/). The engine's workers close after 3 minutes without a calculation, one set kept ready with the modules
  used, and reopen at the first gesture (a field clicked into or changed, a stage made current, the tab back on
  one); the twin's shaders are compiled when Analyze opens, not when the app does. Plan: "1 point" of power, and the
  Assumed TGI chart stays hidden until its sizes come (it showed a missing column meanwhile).

* Plan, second pass (2026-11-11): the animals per group and the end of study leave the form: a click on the
  Analyzable or End of study card sets them, the card says "set", and another click changes them or goes back to the
  computed value. Multiple comparisons moves to Advanced (Holm by default, none needed with one arm; the protocol names
  it). The calendar's instructions show while the pointer is on it, the first three times. When the controls may level
  off below the endpoint and no lower endpoint avoids it (syngeneic lines), the warning says to plan to end on the
  study's day instead of suggesting a lower endpoint (that day is the plan's end in every study: Crux's endpoint day
  is the median over all studies, those whose controls never reach it counted as never, tgi design.rs:1130).
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
  `8e2184ef…b71a6c` with LF line endings, the one PMX Studio records and the git blob; `e2c069fe…d605` for the
  Windows CRLF checkout, as in the seal manifest and its addendum
  `analysis/out_regen/SEAL_MANIFEST_2026-11-11_ADDENDUM_2026-11-13.txt`), in place of the export of 2026-11-06: the curation notes applied and three exclusions of
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
* Overview: the app opens on its stages as cards (Plan, Monitor, Analyze) joined by the files that pass between them, the user's studies (saved plans and followed studies), and Open a study: a plan, a followed study or a table of measurements, dropped anywhere or picked, opens in its stage. `/plan` opens Plan directly
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
