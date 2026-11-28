# PMX Studio: developer notes

The user-facing help is [README.md](README.md). This file is for developers: how the package is built, its
functions, files and formats.

## Architecture

> **Prototype.** This package is the prototype of PMX Studio Preclinical. Its study design and analysis run on
> [Crux](https://github.com/datagrok-ai/crux): `stats.tumor_study_design`, `stats.tgi`, `stats.tumor_twin` and
> `stats.tumor_body_weight` (crate tgi), with every other calculation also a Crux call. The package holds the user interface (U2, `libraries/u2`), the cell lines'
> priors and the therapies' archetypes, and turns Crux's results into Datagrok tables. Its growth and weight
> charts are a band chart of its own (`src/band-chart.ts`, on echarts): the platform's line chart cannot shade the
> area between two columns.

## What the engine computes

De novo design of xenograft and syngeneic efficacy studies: the user declares the cell line, the endpoint
volume and the treated arms (each with its kind of therapy and expected TGI or T/C); the design gives the end
of the study and its uncertainty, the analysable animals per arm for the declared power (with a multiplicity
correction over the arms, Bonferroni by default: each arm against the control by Welch's t test), the animals to implant for them, the power curves and the
predicted volumes of every arm; given the animals per arm, the power there and each arm's smallest detectable
effect. Analyze reads a finished study as Plan pre-specified it (each arm against the control on the per-animal
slope of log volume, with the same correction) and adds the digital twin of every treated animal: what it would
measure untreated, from its volume at randomization, the study's controls and the cell line's model; its paired
effect on the evaluation day; and each control checked against the others. Given the plan saved in Plan, Analyze
also reads the study against it: whether the controls grew as predicted, whether each treated arm's effect was the
one assumed, and whether the study ended when predicted. Given the animals' body weight, Analyze also reads whether
the treatments were tolerated.

The control arm comes from a curated atlas of historical controls mined from open-access supplementary
data with [unpool](https://github.com/gcanudo-barreras/unpool) (12 tumor models: 7 syngeneic lines, 3 cell-line xenografts and 2 PDX), fitted as a reduced
Gompertz model per model; a new study's growth is drawn from the atlas's distribution of labs (`study_prior` in the
priors file).

## Functions

| Function | What it returns |
|---|---|
| `pmxPriors` | the control-growth prior of each cell line: initial growth rate alpha (new-study median and 95% interval), deceleration beta, between-study tau, ESS, non-take rate, papers and animals |
| `pmxPlanStudy(line, endpoint, arms, power, correction, n?, lastDay?, tgiScenarios?, alpha?, preLoss?, postLoss?, lossRule?, confidence?, vRand?, measurementDays?, startWeekday?, humane?, compareSchedules?, fits?, effectAs?)` | eight tables (`arms`: name, therapy, metric, effect, and optionally `therapy_source`, `name` when the form reads the therapy from the name; `effectAs`, TGI or T/C, sets every arm's metric): `scenarios` (per assumed TGI in `tgiScenarios` and measurements per week: the animals needed, scenarios and not predictions; empty without it), `plan` (the study: its end and uncertainty, the animals to analyse, randomize and implant, the per-comparison level, the measurement days, the notes), `armResults` (per arm: expected and realised effect, kill rate, animals it needs alone, power at the study's size, smallest detectable effect at a given n), `powerCurve` (power per arm and n), `curves` (predicted volumes of the control and every arm, with 50% / 90% bands), `duration` (power and the controls past the endpoint and the humane limit, by the study's last day), `schedules` (power and animals by measurements per week; empty with `compareSchedules` false, the Plan view's choice while that section is folded) and `calendar` (what is done each day from randomization) |
| `pmxAnalyzeStudy(study, animal, group, day, volume, bodyWeight?, tumorWeight?, removalReason?, plan?, line, control, endpoint, correction, evaluationDay?, alpha?, weightLimit?, humane?, density?)` | eleven tables: `analysis` (one row: the evaluation day, the control, the control check's coverage, the notes, and against the plan the controls' and the end's verdicts), `comparisons` (per arm: TGI and T/C on the evaluation day, without an interval; the slope test against the control with its interval at `slope_difference_confidence`, the median paired effect and the clear responders), `animals` (slope, time to endpoint, mRECIST), `volumes` (observed median and range per group and day), `twins` (per treated animal and day: measured volume and its untreated twin's median and 90% band), `effects` (per treated animal: its paired effect on the evaluation day and 90% interval) `checks` (per control measurement: its band from the other controls), `versusPlan` (per arm: the plan arm it ran as, the assumed and the observed TGI or T/C, and whether the assumption held) `predicted` (the controls as predicted and as observed), `tolerability` (per group, as weighed and tumor-corrected: the weight-loss grades, the mean nadir, the median day of the worst loss, the recoveries and the animals removed early, with who and why) and `weights` (per group and day, both readings: the mean body weight change from day 0 and its range); the last two are empty without `bodyWeight`. `plan` is the text of a `.pmxplan` file; empty, the plan the study table carries in its `pmx.plan` tag, else the controls and the end are read against the line's atlas |
| `pmxPlanFileHandler(content)` | opens a `.pmxplan` file in Plan (double-click in Files, or a drop) |
| `pmxMonitorStudy(study, animal, group, day, volume, ...)` | six tables for a running study: `status`, `animals` (each animal's last value, its expected next one and why it needs attention), `days` (each group's forecast), `decisions` (on a planned review), `tolerability` and `weights` |
| `pmxStudyFileHandler(content)` | opens a `.pmxstudy` file in Monitor |
| `pmxPlanFilePreview(file)`, `pmxStudyFilePreview(file)` | a `.pmxplan` or `.pmxstudy` selected in Files: what it holds, read from the file alone (a blinded study's groups not named), and a button to open it in its stage |
| `pmxLines()` | the atlas's tumor models, for the functions' `line` choices |
| `pmxMarkStudyRows`, `pmxStudyTableView` | the study table's menu actions (flag, exclude, clear) and its decoration, called from `detectors.js` |
| `pmxStudioApp` | no path: the home (Overview in the stage steps: Plan, Monitor and Analyze as cards joined by the files that pass between them, the user's saved plans, and Open a study: a `.pmxplan`, `.pmxstudy` or table of measurements, dropped anywhere on the page or picked, read and opened in its stage, the stage suggested from the dates alone so a running study stays blind); `/analyze`: the Analyze view (U2): the form generated from `pmxAnalyzeStudy`, the figures, the groups with their result (marked when a group is over the NCI limits for tolerability), the observed growth, each treated animal against its untreated twin (the animals listed by effect, the one selected drawn), Plan vs observed (each arm's effect against the assumed one, the controls against their predicted growth, the end against the predicted day, and "Contribute to the atlas", coming soon: what the study's controls could add to the line's predictions), Tolerability (each group's body weight change from day 0, as weighed or tumor-corrected, the weight-loss grades, the animals removed early and each group against the NCI limits; without a body weight column, how to add one), then the results in words ready to copy, the statistics, the control check and the analysis protocol ready to copy; `/plan`: the Plan view (U2): the declaration form generated from `pmxPlanStudy`, the arms as cards, the results and the study calendar (its measurement and randomization weekdays set on it), the alternatives folded, Explore (the TGIs entered in the form, each a scenario and not a prediction, and the animals each needs), the study protocol ready to copy, and Save plan / Download in the ribbon |

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
| `compare-view.ts` | the Compare view: studies side by side, each control against the atlas |
| `handlers.ts`, `pmx-objects.ts` | animals, groups and figures as Datagrok objects: their handlers and actions |
| `tolerability.ts` | its Tolerability tab |
| `plan-file.ts` | the saved plan (`.pmxplan`): built from Plan's call and outputs, checked when read, saved to Home or downloaded |
| `ui.ts` | what the views share: formats, blocks, notes, the volume axis, a stage's two panes and the recompute on every change (no Run) |
| `stages.ts` | moving between the home and the stages' views |
| `hub-view.ts`, `css/hub.css` | the home and its Open a study dialog |
| `open-study.ts` | what a file brought to the home holds and the stage it opens in |
| `monitor.ts`, `monitor-view.ts` | a running study through Crux `stats.tumor_monitor` and `stats.tumor_allocate`, and the Monitor view |
| `study-file.ts`, `study-meta.ts` | the followed study (`.pmxstudy`) and the study's details (ID, owner, day 0, what each group receives, notes) |
| `study-table.ts`, `outlier-rules.ts`, `detectors.js` | the study's own table view: the marks (Excluded, Flagged), the outlier rules, paste, the context menu |
| `study-header.ts`, `panels.ts` | the study's header and stage steps; a stage's panels, tabs and layout editor |
| `library.ts` | the Overview's list of studies |
| `calendar.ts`, `sheet.ts`, `report.ts`, `arrive.ts` | the study calendar and its `.ics`, the randomization sheet, the report and exports, the ARRIVE checklist |
| `docx.ts`, `ethics-documents.ts` | Word files written by hand (a zip of XML parts); the ethics application documents (study plan after the ARRIVE study plan, EU non-technical summary) |
| `survival.ts`, `km-chart.ts`, `spider-chart.ts`, `synergy.ts`, `explore.ts` | days to the endpoint (Kaplan-Meier, projection, Fine-Gray), each animal's plots, combination synergy, the study as tables |
| `plan-vs-observed.ts`, `plan-checks.ts` | the study against its plan, and the plan's predictive ranges |
| `line-input.ts`, `band-chart.ts`, `gpu.ts` | the tumor model picker, the band chart, the WebGPU device |

## Opening a study

**Open a study** on the Overview takes a `.pmxplan`, a `.pmxstudy` or a table of measurements (`.csv`, `.tsv`,
`.xlsx`), browsed for or dropped anywhere on the page. Plan, Monitor and Analyze have the same dialog as **Open…** in
their header, and take files dropped on their page too: Plan opens a plan, or a table as the growth history of a model
outside the atlas; Monitor a followed study, or a running study's table with its plan; Analyze a finished study's
table with its plan (its mouse PK is added in Analyze once it is open). A file of another stage opens in the stage it
belongs to, and the dialog says so. Once read, the dialog shows the study as a summary (groups, animals, days, the files
it came from), what was read and how, and a warning for anything it could not read; the main button opens it where the
dates say it belongs, the others in the other stages it may open in.

Every dialog opens once: a second click while it is open brings the open one to the front (`ui.ts` `modal`).

Opening PMX Studio leaves Datagrok's Browse tree, toolbox and help as they were. Pin the PMX Studio view to give it
the whole screen: Browse closes while the pinned view is current (about a second after you pin it) and comes back as
it was when you leave it or unpin it. The context panel stays closed in the app's views (F4 opens it). The Overview's toolbox has the app's actions:
Open a study, New plan and Compare saved studies (all your studies with measurements; with fewer than two, it shows
Your studies); the demo data; and this help.

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
  reach it, their share said). One test, one interval: the result and its uncertainty come only from the slopes,
  shown as "Growth rate vs control: grows 4.1%/day slower (95% CI 2.0–6.3%), p = 0.003" (the log-slope difference d
  and its adjusted interval as 100 (1 − e^d) per day; "faster" when d > 0; the level is 1 − α/k with a correction).
  TGI and T/C are the evaluation day's descriptions, without an interval (a one-day interval could contradict the
  test, and covers only the animals' sampling, not how labs differ). With no plan, the controls and the end are read against
  the line's atlas with this study's controls.
- Plan's kill rate makes the declared effect the median, over the labs the atlas draws, of the effect each one shows
  when its own controls reach the endpoint (Crux, decision tumor-plan-checks).
- The package's tests use an example study (`src/tests/example-study.ts`, not in the app), run as Plan designs it:
  each arm with the animals Plan asks for (MC38, TGI 40% and 70%: 26 controls, 26 and 6 treated), one lab of the atlas
  drawn with the design's seed and nothing picked by its outcome, each group with mice of its own. (Crux's "median
  study", the lab whose controls reach the endpoint on the median day, is not median in effect: its 40% arm shows a
  TGI of about 21%, so the example no longer uses it.)

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
  102:1555), at the tumor's when its volume reached **Humane limit, tumor** (1500 mm³ unless the lab set its own in the package's settings, `humaneLimit`), else for another reason. With
  a **Removal reason** column (found dead, moribund, weight, tumor, ulceration, other; on any of the animal's rows)
  the reason given replaces the one read.
- **NCI limits** for a treatment in mice (Plowman et al. 1997, in Teicher (ed.) *Anticancer Drug Development
  Guide*): a group's mean loss at the nadir over 20% as weighed, or more than 10% of its animals dead or removed for
  their weight (`NCI` in `analyze.ts`; Crux `toxic_nadir`, `toxic_deaths`; some labs use 20% deaths, Corbett et al.
  1996). Shown per group, the vehicle included; a descriptive reading of an efficacy study, not a safety assessment
  (no NOAEL or MTD), and no animal is excluded from the efficacy analysis for it.
- Doses held are not read yet: a study table has no column for them.
- The test example's body weight is illustrative (Crux `stats.tumor_body_weight_simulate`: a 20 g mouse, a slow gain,
  a dip with its nadir on day 5 of 0%, 4% and 10% for the vehicle and the 40% and 70% arms, the tumor's mass at
  MC38's density with each animal's own spread, removal at 20% loss), and so are its tumor weights, one per animal
  after its last measurement: the atlas has no body weights yet.

## Analyze: exposure and effect

With NCA's exposure kept for the study (its AUClast per group), Analyze's Results tab has an **Exposure and effect**
panel of its own: each treated group's TGI against its mouse exposure, or against its dose when no exposure was kept.
The panel is hidden while there is nothing to draw. Given the patients' exposure of a drug (measured or projected,
typed in Compare or in the panel), the groups table gains a **× patients' exposure** column: the group's mouse AUC as
a multiple of the patients', its tooltip saying which AUC was used, whether it was corrected for protein binding, and
where the patients' figure came from. Exposures, AUCs, Cmax, concentrations, doses and multiples are shown to 2–3
significant figures (`sig` in `sheet.ts`), with thousands separators.

**Animals by day**: the growth chart's tooltip gives each point's animals ("8 of 10 animals"), and a stretch where
fewer than 3 remain is dashed (a median of one or two animals is no group's). Under the chart, one line names the
groups that lost animals. The full per-day table is in the group's context panel, in the report (when a group lost
animals) and under **Data > Open as table > Animals by day**.

## Compare

Compare (the Overview's **Compare saved studies**, or two demo studies) puts studies with measurements side by side,
in panels laid out as a stage's (**View > Layout…**):

- **The studies**: how many, of which models, their treated groups, how many are placed against the patients'
  exposure, and how many controls grew within the atlas. When the studies are of more than one model, a line says that
  tumour growth inhibition is not directly comparable across models.
- **Tumor growth**: every group's median and range (smallest to largest tumor) by day, coloured by its treatment, each
  study after the first dashed; **Linear**, **Log**, or **Relative V/V0** (each animal over its first volume, for
  studies randomized at different sizes).
- **Dose and effect**: each drug alone at each dose it was given, across the studies (a drug given in two units: two
  curves). A dose is read from the group's details, or its name.
- **Exposure and effect** (the title of Analyze's
  panel too): each treated group's TGI against its multiple of the patients' exposure, the
  patients' range shaded, with the same wording as Analyze. A drug without the patients' exposure has an **Add
  patients' exposure…** button, which saves it with its study.
- **Controls against the atlas**: each study's median control over the atlas's prediction for its model (the 95%
  range Analyze draws in grey, randomized at the study's day-0 median), and a table: day 0, the control's doubling
  time (its animals' mean growth rate of log volume), the atlas range, and whether it is within it or how many times
  faster or slower. A control outside the range gets one plain sentence. Studies record no passage, so there is no
  passage column.
- **Treated groups**: study, model, group, TGI (1 − ΔT/ΔC: the analysis's, or the same from the rows, study-effect.ts growthTgi), day, ×
  patients' exposure, the analysis's result when the study kept it, and where the patients' exposure came from.

**Report (PDF)…** prints the page with the exposure section.

## Toolbox

Plan, Monitor, Analyze and Compare have the platform's toolbox on the left. It holds the same items as the view's menu,
grouped as **Study** (Open…, Study details…, Save…; Monitor's Propose groups… and the randomization sheet),
**Export**, **Documents** (ARRIVE, the CTD summary, the ethics documents), **View** (the layout, the tables to open),
and **Help** (the demo, the example study, PMX Studio's help). Analyze's Export also has the study file (`.pmxstudy`)
and the calendar (`.ics`).

## Right-click and the context panel

Animals, groups and figures are Datagrok objects (`handlers.ts`: a `DG.ObjectHandler` for each of `PmxAnimal`,
`PmxGroup` and `PmxFigure`, in `pmx-objects.ts`). Their actions are registered once, with `registerParamFunc`, so a
right click on one and the context panel's **Actions** pane offer the same ones, and an action shows only where it
applies:

- **An animal** (Analyze's twin list, waterfall bars, Monitor's lists): Flag this animal…,
  Exclude from the analysis… or Include in the analysis, Clear the flag, Select the animal's rows, Copy the animal's
  name.
- **A group** (its chip in the study header, its name in the groups table): Change color…, Flag every animal of the
  group…, Exclude the group from the analysis…, Clear its marks, Select the group's rows, Copy the group's name.
- **A figure** (a chart, or a stage's figure cards): Save as PNG, Copy as image; the cards: Pin these figures to
  compare, or Unpin.

The details a click shows (the animal's readings, the group's) are the object's context panel, with its main action
as a button. A chart's bar or a legend item cannot be bound to an object: its right click opens a menu of the same
actions (a legend item adds Change color…). The colour swatches and the calendar's days keep their small menus (a
swatch's Change color… is also the group's).

## Where plans and studies are saved

Plans and studies are saved in PMX Studio itself, as rows of its domain tables (`databases/pmx/schema.json`,
`src/store.ts`): Save (or Ctrl+S) updates the row, the Overview lists them (Delete, Show deleted), Browse lists Plans
and Studies under the app, and `/monitor?study=<id>` or `/plan?plan=<id>` opens one. Plans and studies saved as files
before are brought in once with **Bring your files in** on the Overview; the files are left where they are, and
.pmxplan / .pmxstudy files still open, drop and export as before. A study is yours until you share it.

- Share a plan or study from its row on the Overview (the share icon): the people or groups you pick can open it
  and read it. Only you can change it; to let someone else edit, export the file or ask the owner.
- A study shared with you to read is saved as a copy of your own (**Save as my copy**), marked as a copy of the
  owner's.
- Blinding hides the allocation but doesn't lock it: anyone with read access can see it outside the app, in
  Datagrok's generic table views. Earlier versions are not kept: plans and studies have no change history.
- The row (schema 1.5.0): `document` is the .pmxplan / .pmxstudy document as `json` (jsonb gives its keys back
  reordered, so nothing compares its text: a plan is known by `planKey`, its name and inputs with sorted keys, and
  kept tables by `sameAsked`); `owner` is the user the row belongs to and `copy_owner` the original row's user (both
  `user`, listed by name through query captions); `copy_file_owner` is the owner text of someone else's file brought
  in as a copy; `run_by` is Study details' free-text Run by; `start` is day 0 (`datetime`, midnight UTC: EMS has no
  date type) and `end_day` / `last_day` / `next_day` whole days (`int`, rounded on write). There is no project:
  EMS links no row to a Space (a project in an older document is read as a tag).
- Schema changes go with a migration: 1.4.1 → 1.5.0 changes column types, which the platform does not do on its own,
  so `databases/pmx/0001_schema_1_5_0.sql` (and `down/`) does it, on the pattern of Chem's mpo and Grit's. Its
  `ems-migration` header needs the snapshot hashes from the platform's schema snapshot tool before the package is
  deployed where 1.4.1 rows exist. The app also reads a row as 1.4.1 wrote it (text document and day 0).
- A lab that needs several people editing one study should switch the `domainStorage` setting (**Save in
  Datagrok**) off and use shared files. Export (download) the plans and studies first: with it off, the Overview
  lists files only, and says how many are kept in Datagrok storage.

## Files

- `files/engine_priors_v3.json`: the priors, exported offline from the curated atlas, with the atlas hash and the
  Crux commit they came from; per line, `twin` is the digital twin's own fit and its leave-one-control-out
  validation. v3 keeps v2's MC38, CT26 and B16F10 verbatim and adds the syngeneic 4T1, LLC, Panc02, Hepa1-6 and EMT6,
  the cell-line xenografts MDA-MB-231 and A549 (R7 freeze: SKOV3 out), and PDX-NSCLC and PDX-PDAC (`provenance.v3`: each line's
  tier, and which fields are its own or pooled).
- `files/tumour_density_v0.json`: the tumor's density per line and its spreads (between labs, groups and animals),
  exported offline from the atlas papers' excised tumors.
- `files/SOURCES.md`: the papers the priors and the density were estimated from, with DOI and licence. The per-paper table in the
  priors (`studies`) uses anonymous ids.
- `files/example/mc38-ld4172-anti-pd1.csv`: the example study the empty Monitor and Analyze offer to download (see
  below), its source, licence and changes in `files/example/README.md`.
- `files/demo/`: demo data to try every stage (see "Demo data" below and `files/demo/README.md`).
- `crux-modules/`: Crux's browser modules, the crux-core commit they come from in `crux-modules/VERSION` (see below).

## Example study

The empty Monitor and Analyze offer one published study to try: MC38 tumors, vehicle, LD4172 (a RIPK1 degrader),
anti-PD-1 and the two combined, 8 mice per group, tumor volume and body weight on days 0 to 12 from randomization.
From Yu X. et al., "Development of a RIPK1 degrader to enhance antitumor immunity", Nature Communications 15 (2024),
https://doi.org/10.1038/s41467-024-55006-2, Supplementary Fig. 4A–B (Source Data file), licensed CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/). Changed by us: reshaped to one row per animal and day; day
re-expressed as days from randomization (article day 7 = day 0); animal labels assigned by column order; group names
completed with dose, route and schedule from the article's Methods (µg written "ug"); each weight paired with the
tumor volume in the same column position of its group (the article gives no animal identifiers). Volumes and weights
unchanged. Its control group is also among the MC38 studies the growth model was built from (`files/SOURCES.md`).

Get it from the empty Monitor or Analyze, or open the file from the Overview's **Open a study**. In Monitor it can be
opened with a plan of 3 arms (the groups are matched to the arms by order).

What it shows: in Analyze the median control stays below the 1000 mm³ endpoint, so the groups are compared on the
last day, day 12, and the app says so; the combination slows growth (TGI 79%, p = 0.013) while neither drug alone has
a clear effect (TGI 17% and 18%); the body weights stay within the NCI limits (Monitor's Welfare, Analyze's Tolerability).

## Demo data

`files/demo/` holds a small set of studies to try every stage: a plan, a running study with planted welfare cases, a
finished study with removals, mouse PK for the exposure, a second study to compare, a model outside the atlas, and
other table layouts. They are one click away on the Overview: each stage's card has **Try with demo data** at its foot
(Plan opens the CT26 trametinib + anti-PD-1 plan, Monitor the same study at day 9, Analyze the finished study with its
mouse PK when the NCA package is installed), and while **Your studies** is empty, **Or compare two demo studies** saves
the CT26 and MC38 trametinib studies as yours and opens them in Compare. The same actions are in the Overview's
toolbox. On a server the files are also in **Browse > Files > App Data > Pmxstudio > demo**; open one from **Open a
study** or a stage's **Open…**, or drag it onto the app. `files/demo/README.md` says what to open for each stage. All are
synthetic (seeded models, not real animals or real PK) except the published MC38 study, which is `files/example/`.
Exposure needs the NCA package, whose own demo data is in **App Data > Nca**.

## Crux modules

One Crux module, `pmx`: the functions this package calls (crate tgi's `stats.tgi` and `stats.tumor_*`, synergy,
`stats.contrast`, `stats.posthoc`, `stats.dist_eval`, `survival.fine_gray`) and the other functions of the Crux
features they belong to (tgi, synergy, LS-means, post-hoc tests, scipy.stats, competing risks), which Crux's
`nonclinical` and `stats` modules own and `pmx` carries, so the package ships 5.6 MB of wasm instead of their 21 MB. Built from crux-core
(`crux-modules/VERSION` names the commit) and kept in git, as Chem keeps `vendor/crux-sketch`, so the package
builds and publishes without Rust. To update it:

```bash
cd ~/crux-core
python3 tools/parity/wasm-modules/modules.py build --bindgen web --out ~/crux-modules-web --modules pmx
# copy manifest.json and pmx/crux_wasm{.js,_bg.wasm} into crux-modules/; the loader (crux-modules.js and the files
# it imports) is crates/crux-wasm/js/*.mjs renamed to .js, their imports of each other rewritten to match (Datagrok
# serves .mjs as text/plain): copy it again when it changes
```

Needs Rust 1.95 with the `wasm32-unknown-unknown` target, `wasm-bindgen-cli` 0.2.129, binaryen 133 and Node.
