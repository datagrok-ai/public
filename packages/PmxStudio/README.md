# PMX Studio

> **Prototype.** This package is the prototype of PMX Studio Preclinical. Its study design and analysis run on
> [Crux](https://github.com/datagrok-ai/crux): `stats.tumor_study_design`, `stats.tgi`, `stats.tumor_twin` and
> `stats.tumor_body_weight` (crate tgi), with every other calculation also a Crux call. The package holds the user interface (U2, `libraries/u2`), the cell lines'
> priors and the therapies' archetypes, and turns Crux's results into Datagrok tables. Its growth and weight
> charts are a band chart of its own (`src/band-chart.ts`, on echarts): the platform's line chart cannot shade the
> area between two columns.

De novo design of xenograft and syngeneic efficacy studies: the user declares the cell line, the endpoint
volume and the treated arms (each with its kind of therapy and expected TGI or T/C); the design gives the end
of the study and its uncertainty, the analysable animals per arm for the declared power (with a multiplicity
correction over the arms, Dunnett by default), the animals to implant for them, the power curves and the
predicted volumes of every arm; given the animals per arm, the power there and each arm's smallest detectable
effect. Analyze reads a finished study as Plan pre-specified it (each arm against the control on the per-animal
slope of log volume, with the same correction) and adds the digital twin of every treated animal: what it would
measure untreated, from its volume at randomization, the study's controls and the cell line's model; its paired
effect on the evaluation day; and each control checked against the others. Given the plan saved in Plan, Analyze
also reads the study against it: whether the controls grew as predicted, whether each treated arm's effect was the
one assumed, and whether the study ended when predicted. Given the animals' body weight, Analyze also reads whether
the treatments were tolerated.

The control arm comes from a curated atlas of historical controls mined from open-access supplementary
data with [unpool](https://github.com/gcanudo-barreras/unpool) (MC38, CT26, B16F10), fitted as a reduced
Gompertz model per cell line; the MAP prior of a new study's growth is Crux `stats.map_prior`.

## Functions

| Function | What it returns |
|---|---|
| `pmxPriors` | the control-growth prior of each cell line: initial growth rate alpha (new-study median and 95% interval), deceleration beta, between-study tau, ESS, non-take rate, papers and animals |
| `pmxMapPrior(line)` | the MAP prior of a new study's log alpha, recomputed in the browser by Crux `stats.map_prior` |
| `pmxPlanStudy(line, endpoint, arms, power, correction, n?, lastDay?, alpha?, preLoss?, postLoss?, lossRule?, confidence?, vRand?, measurementDays?, startWeekday?, humane?, compareSchedules?)` | seven tables: `plan` (the study: its end and uncertainty, the animals to analyse, randomize and implant, the per-comparison level, the measurement days, the notes), `armResults` (per arm: expected and realised effect, kill rate, animals it needs alone, power at the study's size, smallest detectable effect at a given n), `powerCurve` (power per arm and n), `curves` (predicted volumes of the control and every arm, with 50% / 90% bands), `duration` (power and the controls past the endpoint and the humane limit, by the study's last day), `schedules` (power and animals by measurements per week; empty with `compareSchedules` false, the Plan view's choice while that section is folded) and `calendar` (what is done each day from randomization) |
| `pmxAnalyzeStudy(study, animal, group, day, volume, bodyWeight?, tumorWeight?, removalReason?, plan?, line, control, endpoint, correction, evaluationDay?, alpha?, weightLimit?, humane?, density?)` | eleven tables: `analysis` (one row: the evaluation day, the control, the control check's coverage, the notes, and against the plan the controls' and the end's verdicts), `comparisons` (per arm: TGI with its 95% interval, T/C, the slope test against the control, the median paired effect and the clear responders), `animals` (slope, time to endpoint, mRECIST), `volumes` (observed median and range per group and day), `twins` (per treated animal and day: measured volume and its untreated twin's median and 90% band), `effects` (per treated animal: its paired effect on the evaluation day and 90% interval) `checks` (per control measurement: its band from the other controls), `versusPlan` (per arm: the plan arm it ran as, the assumed and the observed TGI or T/C with its 95% interval, and whether the assumption held) `predicted` (the controls as predicted and as observed), `tolerability` (per group, as weighed and tumor-corrected: the weight-loss grades, the mean nadir, the median day of the worst loss, the recoveries and the animals removed early, with who and why) and `weights` (per group and day, both readings: the mean body weight change from day 0 and its range); the last two are empty without `bodyWeight`. `plan` is the text of a `.pmxplan` file; empty, the plan the study table carries in its `pmx.plan` tag, else the controls and the end are read against the line's atlas |
| `pmxPlanFileHandler(content)` | opens a `.pmxplan` file in Plan (double-click in Files, or a drop) |
| `pmxStudioApp` | `/analyze`: the Analyze view (U2): the form generated from `pmxAnalyzeStudy`, the figures, the groups with their result (marked when a group is over the NCI limits for tolerability), the observed growth, each treated animal against its untreated twin (the animals listed by effect, the one selected drawn), Plan vs observed (each arm's effect against the assumed one, the controls against their predicted growth, the end against the predicted day, and "Contribute to the atlas", coming soon: what the study's controls could add to the line's predictions), Tolerability (each group's body weight change from day 0, as weighed or tumor-corrected, the weight-loss grades, the animals removed early and each group against the NCI limits; without a body weight column, how to add one), then the results in words ready to copy, the statistics, the control check and the analysis protocol ready to copy; an example study (which carries its plan) analyzed in one click; any other path: the Plan view (U2): the declaration form generated from `pmxPlanStudy`, the arms as cards, the results and the study calendar (its measurement and randomization weekdays set on it), the alternatives folded, the study protocol ready to copy, and Save plan / Download in the ribbon |

## Sources (`src/`)

| File | Role |
|---|---|
| `engine.ts` | translates a declaration into Crux `stats.tumor_study_design` calls (one per arm, the same seed) and combines them: the study's n is the largest arm's |
| `archetypes.ts` | the treated arm per kind of therapy (onset, between-animal spread of the kill rate), estimated from the atlas's treated arms |
| `priors.ts` | the priors file and the MAP prior through Crux |
| `crux.ts` | Crux's wasm modules in the browser (leaf ABI) |
| `frames.ts` | the design as Datagrok tables (units and reference lines included) |
| `plan-view.ts`, `arms-input.ts` | the Plan view and its arms and multiple-comparisons editors |
| `analyze.ts` | translates a study into Crux `stats.tgi`, `stats.contrast`, `stats.tumor_twin` and `stats.tumor_body_weight` calls and combines them |
| `analyze-view.ts` | the Analyze view |
| `tolerability.ts` | its Tolerability tab |
| `plan-file.ts` | the saved plan (`.pmxplan`): built from Plan's call and outputs, checked when read, saved to Home or downloaded |
| `ui.ts` | what the two views share: formats, blocks, notes, the volume axis, the panes and the stage strip |
| `stages.ts` | moving between the stages' views |

## Saved plans

A plan is a `.pmxplan` file, JSON: `format` `pmx-plan`, `version` 1, `name`, `created`, `author`, `inputs` (the
Plan form's values, to reopen it) and `predicted` (what Plan showed, frozen: line, endpoint, volume at
randomization, end of study and endpoint day with 8 in 10 studies' range, each arm with its metric, effect and
animals, and `checks`: what a study run as planned shows, from Crux `stats.tumor_plan_checks` with the animals each
group has, as 95% ranges over simulated labs: the day the median control reaches the endpoint, the controls' mean
growth, each arm's TGI and T/C, and a study's median control per day). Analyze compares a study with `predicted`,
never with a plan recomputed on a later atlas. The checks are computed when the plan is saved or downloaded, not on
every replan.

- **Save plan** writes it to `<user>:Home/PMX Studio/plans/`; from Files it is shared like any file, and opens in
  Plan on a double-click.
- **Download** gives the same file to send to someone outside Datagrok; dropped on Datagrok, it opens in Plan.
- In Analyze, the **Plan** field takes a path in Datagrok's files, a file from the computer or a drop, and lists the
  user's saved plans; picking one writes its pre-specified analysis (cell line, endpoint, correction, α) into the
  form. A treated group is matched to the plan's arm by name, else by order when the counts agree (said in a note).
- Each treated group's result comes first, in both tabs, from the pre-specified test on the per-animal slopes of
  log volume: "Slows growth" when its interval lies below 0, "Grows faster" above, else "No clear effect" (4.9%
  false alarms in 2000 simulated MC38 studies with no effect, one treated arm, 8 animals per group). For
  a group with no clear effect the tab says what the study could have detected: with a plan, the power the plan
  gave the arm with the animals it had and the animals the plan needed (a no-effect result from an underpowered
  study is not evidence of inactivity); with no plan, the smallest TGI those animals detect in 80 of 100 labs
  (computed only for such groups: a search of a few seconds).
- Each verdict against the plan compares an observed figure with the range of a study run as planned, so a study
  that went as planned is flagged about one time in twenty: an arm's effect is "As planned" when its TGI (or T/C)
  falls in its range and the test tells the group from the control, "Inconclusive" when it falls there but the test
  does not, else stronger or weaker; the controls are "As predicted" when their mean growth (the slope of log volume the test
  compares, shown as a daily percent) falls in its range, else faster or slower; the end is "As predicted" when the
  median control reaches the endpoint within its range (with no top when more than 2.5% of such studies never
  reach it, their share said). A study's own 95% interval is shown but does not decide:
  it covers its animals' sampling, not how labs differ (in 300 on-plan MC38 studies with 8 animals per group, the
  assumed TGI fell outside the observed interval in 26%). With no plan, the controls and the end are read against
  the line's atlas with this study's controls.
- Plan's kill rate makes the declared effect the median, over the labs the atlas draws, of the effect each one shows
  when its own controls reach the endpoint (Crux, decision tumor-plan-checks).
- The example study is run as Plan designs it: each arm with the animals Plan asks for (MC38, TGI 40% and 70%:
  26 controls, 26 and 6 treated), one lab of the atlas drawn with the design's seed and nothing picked by its
  outcome, each group with mice of its own. (Crux's "median study", the lab whose controls reach the endpoint on
  the median day, is not median in effect: its 40% arm shows a TGI of about 21%, so the example no longer uses it.)

## Tolerability

With a **Body weight** column (grams, weighed with the tumor, on the rows of the measurements), Crux
`stats.tumor_body_weight` reads each animal's weight from its baseline (its last weighing on or before day 0):

- **Tumor-corrected**: each weighing less the tumor's mass (a 1500 mm³ tumor is over 1 g of a 20 g mouse and hides
  a loss of the animal), its caliper volume times the line's typical density: the median over the atlas papers that
  weighed their excised tumors (`files/tumour_density_v0.json`; MC38 0.75, CT26 0.92, B16F10 1.07 mg per mm³; about
  ±35% between labs), or **Tumor density** (Advanced) when given. With a **Tumor weight** column (grams, the excised
  tumor on the row of its last measurement), the excised tumors of 200 mm³ or more set the study's density and each
  group's (a normal hierarchy on log mass / volume from the line to the study to its groups, Crux's `lab_sd`,
  `group_sd` and `within_sd` from the same file), and on that row the tumor's own weight is taken off. The tab's
  switch shows the chart, the figures and the grades either way; its note says which density it used.
- **Grades G0–G3** of each animal's largest loss: below 5%, 5–10%, 10–20%, 20% or more, the CTCAE v5 weight-loss
  grades as preclinical oncology adapts them (in mice G3 is the usual humane endpoint); fixed, so that studies
  compare. **Recovered**: a loss of 5% or more back within 5% of day 0 at the last weighing. **Mean nadir**: a
  group's lowest mean change over the days at least half of it was weighed, and a treated group's difference from
  the vehicle's in percentage points (a tumor-bearing vehicle can lose weight itself: CT26 is cachectic).
- **Removed early**: weighed for the last time before the rest of its group; at the weight endpoint when its
  measured loss there reached **Humane limit, weight loss** (Advanced, 20%: Workman et al. 2010, Br J Cancer
  102:1555), at the tumor's when its volume reached **Humane limit, tumor** (2000 mm³), else for another reason. With
  a **Removal reason** column (found dead, moribund, weight, tumor, ulceration, other; on any of the animal's rows)
  the reason given replaces the one read.
- **NCI limits** for a treatment in mice (Plowman et al. 1997, in Teicher (ed.) *Anticancer Drug Development
  Guide*): a group's mean loss at the nadir over 20% as weighed, or more than 10% of its animals dead or removed for
  their weight (`NCI` in `analyze.ts`; Crux `toxic_nadir`, `toxic_deaths`; some labs use 20% deaths, Corbett et al.
  1996). Shown per group, the vehicle included; a descriptive reading of an efficacy study, not a safety assessment
  (no NOAEL or MTD), and no animal is excluded from the efficacy analysis for it.
- Doses held are not read yet: a study table has no column for them.
- The example's body weight is illustrative (Crux `stats.tumor_body_weight_simulate`: a 20 g mouse, a slow gain,
  a dip with its nadir on day 5 of 0%, 4% and 10% for the vehicle and the 40% and 70% arms, the tumor's mass at
  MC38's density with each animal's own spread, removal at 20% loss), and so are its tumor weights, one per animal
  after its last measurement: the atlas has no body weights yet.

## Files

- `files/engine_priors_v2.json`: the priors, exported offline from the curated atlas, with the atlas hash and the
  Crux commit they came from; per line, `twin` is the digital twin's own fit and its leave-one-control-out
  validation.
- `files/tumour_density_v0.json`: the tumor's density per line and its spreads (between labs, groups and animals),
  exported offline from the atlas papers' excised tumors.
- `files/SOURCES.md`: the papers the priors and the density were estimated from, with DOI and licence. The per-paper table in the
  priors (`studies`) uses anonymous ids.
- `crux-modules/`: Crux's browser modules, the crux-core commit they come from in `crux-modules/VERSION` (see below).

## Crux modules

The two Crux modules the package's functions are routed to (`nonclinical`: crate tgi's `stats.tumor_*`; `stats`:
`stats.contrast`, `stats.dist_eval`, `stats.map_prior`), built from crux-core main (`crux-modules/VERSION` names the
commit) and kept in git, as Chem keeps `vendor/crux-sketch`, so the package builds and publishes without Rust. To
update them:

```bash
cd ~/crux-core
python3 tools/parity/wasm-modules/modules.py build --bindgen web --out ~/crux-modules-web --modules nonclinical,stats
# copy manifest.json, the loader's *.mjs renamed to .js with their imports of each other rewritten to match (Datagrok serves .mjs as text/plain), and <module>/crux_wasm{.js,_bg.wasm} of nonclinical and stats into crux-modules/
```

Needs Rust 1.95 with the `wasm32-unknown-unknown` target, `wasm-bindgen-cli` 0.2.129, binaryen 133 and Node.
