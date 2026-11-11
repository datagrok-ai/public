/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PlanFrames, declarationOf, mapPriorFrame, planFrames, priorsFrame, tgisOf} from './frames';
import {openPlanFile, planView} from './plan-view';
import {AnalysisFrames, analysisDeclarationOf, analysisFrames, warmTwin} from './analyze';
import {analyzeView} from './analyze-view';
import {MonitorFrames, monitorFrames} from './monitor';
import {monitorView, openStudyFile} from './monitor-view';
import {hubView} from './hub-view';
import {quietContextPanel} from './stages';
import {decorateStudyView, markRows} from './study-table';
import {tumorModels} from './priors';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.app({name: 'PMX Studio',
    description: 'Plans, follows and analyzes xenograft and syngeneic efficacy studies.'})
  static async pmxStudioApp(
    // the app's URL path, which the platform passes: /plan opens Plan, /analyze Analyze, anything else the home
    @grok.decorators.param({options: {metaUrl: true, optional: true}}) path?: string): Promise<DG.ViewBase> {
    quietContextPanel();
    const stage = path?.replace(/^\//, '') ?? '';
    // opened in the background, ready for the first calculation that wants it; the twin's shaders compiled in Crux's
    // worker (a failure only leaves the first twin slower)
    void warmTwin().catch(() => null);
    return stage.startsWith('analyze') ? analyzeView() : stage.startsWith('monitor') ? monitorView() :
      stage.startsWith('plan') ? planView() : hubView();
  }


  @grok.decorators.fileHandler({ext: 'pmxstudy', description: 'Opens a PMX Studio study in Monitor.'})
  static async pmxStudyFileHandler(content: string): Promise<DG.DataFrame[]> {
    quietContextPanel();
    return openStudyFile(content);
  }

  @grok.decorators.func({
    description: 'The tumor models in the atlas of historical controls (syngeneic lines and PDX): the choices of every ' +
      '`line` input',
    outputs: [{name: 'result', type: 'list<string>'}],
  })
  static async pmxLines(): Promise<string[]> {
    return Object.keys(await tumorModels());
  }

  @grok.decorators.func({
    name: 'pmxMonitorStudy',
    description: 'Follows a running xenograft or syngeneic efficacy study with Crux stats.tumor_monitor, from every measurement so far, the cell line’s growth model from the atlas and the plan: for each animal the range of its next measurement, its chance of passing the humane limit and whether its last value was unexpected; for each group its median day by day to the end; the day the study is expected to end; and, on the reviews fixed at the start, each treated group’s chance that the final test finds it slower than the control.',
    outputs: [
      {name: 'status', type: 'dataframe', options: {description: 'The study today: one row.'}},
      {name: 'animals', type: 'dataframe', options: {description: 'Per animal: the range of its next measurement and why it needs attention.'}},
      {name: 'days', type: 'dataframe', options: {description: 'Per group and day, measured and to come: its median volume and range.'}},
      {name: 'tolerability', type: 'dataframe', options: {description: 'Per group: its body weight, as Tolerability reads it (empty without body weights).'}},
      {name: 'weights', type: 'dataframe', options: {description: 'Per group and day: its mean body weight change (empty without body weights).'}},
      {name: 'decisions', type: 'dataframe', options: {description: 'Per treated group, on a planned review: the chance of success (empty between reviews).'}},
    ],
  })
  static async pmxMonitorStudy(
    @grok.decorators.param({options: {category: 'Study', caption: 'Study', nullable: false, description: 'The measurements so far: one row per animal and day, day 0 the randomization.'}}) study: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Animal', description: 'The animal (tumor) identifier, unique in the study.'}}) animal: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Group', description: 'The treatment group.'}}) group: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Day', type: 'numerical', description: 'Days from randomization.'}}) day: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Volume', type: 'numerical', description: 'Tumor volume (mm³).'}}) volume: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Body weight', type: 'numerical', optional: true, nullable: true, description: 'Body weight (g), weighed with the tumor (empty: no weight checks).'}}) bodyWeight: DG.Column | null,
    // @ts-expect-error (the platform's `visible` rule, as in pmxPlanStudy)
    @grok.decorators.param({type: 'string', options: {category: 'Study', caption: 'Plan', optional: true, nullable: true, visible: 'false', description: 'The plan saved in Plan (its .pmxplan file’s contents).'}}) plan: string | null,
    @grok.decorators.param({options: {category: 'Study', caption: 'Planned reviews', optional: true, nullable: true, description: 'The study days fixed at the start to review each treated group, comma-separated, e.g. 7, 10.'}}) looks: string | null,
    @grok.decorators.param({options: {category: 'Monitoring', caption: 'Cell line', choices: 'PmxStudio:pmxLines()', initialValue: 'MC38', description: 'The tumor model: its growth model from the atlas of historical controls makes the forecasts.'}}) line: string,
    @grok.decorators.param({options: {category: 'Monitoring', caption: 'Control group', optional: true, nullable: true, description: 'The control (vehicle) group (empty: the one named vehicle, control, PBS, saline or untreated, else the first).'}}) control: string,
    @grok.decorators.param({options: {category: 'Monitoring', caption: 'Endpoint volume', units: 'mm³', min: '500', max: '3000', initialValue: '1000', description: 'The study ends when the median control reaches it.'}}) endpoint: number,
    @grok.decorators.param({options: {category: 'Monitoring', caption: 'Humane limit', units: 'mm³', min: '500', max: '4000', initialValue: '1500', description: 'Your protocol’s limit: an animal near it is listed. Workman 2010 advises a mean diameter of 1.5 cm at most in therapy studies, about 1,700 mm³. Your lab’s limit is set once, in the package’s settings.'}}) humane: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Multiple comparisons', choices: ['holm', 'bonferroni', 'none'], initialValue: 'holm', description: 'Correction over the groups compared with the control in the final test.'}}) correction: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Evaluation day', optional: true, nullable: true, description: 'Leave empty for the day the median control reaches the endpoint volume.'}}) evaluationDay?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Significance (α)', optional: true, min: '0.001', max: '0.2', initialValue: '0.05', description: 'Two-sided significance level of the final test.'}}) alpha?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Humane limit, weight loss', optional: true, units: '%', min: '5', max: '50', initialValue: '20', description: 'An animal is listed when it nears this loss of weight, and is over it at this one (your institution’s limit).'}}) weightLimit?: number,
  ): Promise<MonitorFrames> {
    return monitorFrames({study, animal, group, day, volume, bodyWeight, weightLimit, plan, looks, line, control,
      endpoint, humane, correction, evaluationDay, alpha});
  }

  // the extension as a literal (the decorators' generator reads no constants): plan-file.ts's PLAN_EXT
  @grok.decorators.fileHandler({ext: 'pmxplan', description: 'Opens a PMX Studio plan in Plan.'})
  static async pmxPlanFileHandler(content: string): Promise<DG.DataFrame[]> {
    quietContextPanel();
    return openPlanFile(content);
  }

  @grok.decorators.func({
    name: 'pmxMarkStudyRows',
    description: 'Flags rows of a study table to look at, excludes them from the analysis (the reason asked in a dialog, required to exclude) or clears their marks; PMX Studio\'s items in the grid\'s context menu call it.',
  })
  static pmxMarkStudyRows(
    @grok.decorators.param({options: {description: 'The study table.'}}) table: DG.DataFrame,
    @grok.decorators.param({options: {description: 'The rows, comma-separated indexes.'}}) rows: string,
    @grok.decorators.param({options: {choices: ['flag', 'exclude', 'clear'], description: 'What to do with them.'}}) action: string,
    @grok.decorators.param({options: {description: 'The rows in words, for the dialog and the study\'s log.'}}) what: string): void {
    markRows(table, rows.split(',').map(Number), action, what);
  }

  @grok.decorators.func({
    name: 'pmxStudyTableView',
    description: 'Shows in a study\'s table view what PMX Studio\'s outlier rules find, the marks and their reasons, and PMX Studio\'s menu.',
  })
  static pmxStudyTableView(
    @grok.decorators.param({type: 'view', options: {description: 'A table view of a study.'}}) view: DG.View): void {
    if (view instanceof DG.TableView)
      decorateStudyView(view);
  }

  @grok.decorators.func({
    name: 'pmxPriors',
    description: 'Control-growth priors of the PMX Studio atlas (reduced Gompertz) per syngeneic cell line: study-level initial growth rate alpha (median and 95% interval of a new study), deceleration beta, between-study tau, effective sample size, non-take rate, papers and animals behind each prior.',
  })
  static async pmxPriors(): Promise<DG.DataFrame> {
    return priorsFrame();
  }

  @grok.decorators.func({
    name: 'pmxMapPrior',
    description: 'Recomputes the MAP prior of a new study initial growth rate (log alpha) for a cell line with Crux stats.map_prior from the atlas per-paper estimates: one row per mixture component.',
  })
  static async pmxMapPrior(
    @grok.decorators.param({options: {choices: 'PmxStudio:pmxLines()', description: 'Tumor model (syngeneic cell line or PDX).'}}) line: string): Promise<DG.DataFrame> {
    return mapPriorFrame(line);
  }

  @grok.decorators.func({
    name: 'pmxPlanStudy',
    description: 'Designs a xenograft or syngeneic efficacy study de novo with Crux stats.tumor_study_design: from the cell line, the endpoint volume and the treated arms (each with its kind of therapy and expected TGI or T/C), predicts the end of the study (when the median control reaches the endpoint) and its uncertainty, each arm’s kill rate, the analysable animals per arm for the declared assurance (power over the historical-control prior) with the multiplicity correction chosen (Holm by default; each arm against the control by Welch’s t test), the animals to implant for them, the power curves and the predicted control and treated arms. Given a number of animals per arm, the power there and each arm’s minimum detectable effect. Given an end of study, the study is evaluated that day; the power by end of study shows what ending earlier or later costs, and the power by measurement schedule what measuring less often costs.',
    outputs: [
      {name: 'plan', type: 'dataframe', options: {description: 'The study: one row.'}},
      {name: 'armResults', type: 'dataframe', options: {description: 'The treated arms: one row each.'}},
      {name: 'powerCurve', type: 'dataframe', options: {description: 'Assurance (power averaged over the labs of the atlas) per arm and number of analysable animals per arm.'}},
      {name: 'curves', type: 'dataframe', options: {description: 'Predicted volumes of the control and each treated arm, day by day from randomization: median and 50% / 90% bands.'}},
      {name: 'duration', type: 'dataframe', options: {description: 'Each arm’s power by the day the study ends, with the share of controls past the endpoint and over the humane limit.'}},
      {name: 'schedules', type: 'dataframe', options: {description: 'Each arm’s power and the analyzable animals the study needs, by measurements per week (empty when not asked for).'}},
      {name: 'calendar', type: 'dataframe', options: {description: 'The study day by day from randomization: what is done each day.'}},
      {name: 'scenarios', type: 'dataframe', options: {description: 'Per assumed TGI and measurements per week: the animals the study needs with every arm at that TGI (empty when no TGIs are assumed). Scenarios, not predictions.'}},
    ],
  })
  static async pmxPlanStudy(
    @grok.decorators.param({options: {category: 'Tumor model', caption: 'Cell line', nullable: false, choices: 'PmxStudio:pmxLines()', initialValue: 'MC38', description: 'Tumor model (syngeneic cell line or PDX): its growth prior comes from the atlas of historical controls.'}}) line: string,
    @grok.decorators.param({options: {category: 'Tumor model', caption: 'Endpoint volume', nullable: false, units: 'mm³', min: '500', max: '3000', initialValue: '1000', description: 'The study ends when the median control reaches it.'}}) endpoint: number,
    @grok.decorators.param({options: {category: 'Treatment arms', caption: 'Arms', optional: true, description: 'The treated arms, each compared with the control: name, therapy (small_molecule, biologic, immunotherapy, immune_agonist, cell_therapy, vaccine_peptide, oligo_rna), metric (TGI or TC) and the expected effect in percent. Empty: one small-molecule arm at TGI 60%.'}}) arms: DG.DataFrame,
    @grok.decorators.param({options: {category: 'Design', caption: 'Target power', nullable: false, units: '%', min: '50', max: '99', initialValue: '80', description: 'Chance of detecting each arm’s effect.'}}) power: number,
    @grok.decorators.param({options: {category: 'Design', caption: 'Multiple comparisons', nullable: false, choices: ['holm', 'bonferroni', 'none'], initialValue: 'holm', description: 'Correction over the arms compared with the control, each by Welch’s t test: Holm (recommended), Bonferroni, or none (each comparison at alpha, as in an exploratory study).'}}) correction: string,
    @grok.decorators.param({type: 'int', options: {category: 'Design', caption: 'Animals per arm', optional: true, nullable: true, min: '2', max: '100', description: 'Empty for the number the target power needs, or set it to see the power and smallest detectable effect with that many.'}}) n?: number,
    @grok.decorators.param({type: 'int', options: {category: 'Design', caption: 'End of study', optional: true, nullable: true, units: 'day', min: '5', max: '120', description: 'Leave empty to end when the median control reaches the endpoint.'}}) lastDay?: number,
    // @ts-expect-error (the platform's `visible` rule, as above)
    @grok.decorators.param({type: 'string', options: {category: 'Explore', caption: 'Assumed TGIs', optional: true, nullable: true, units: '%', visible: 'false', description: 'Up to 8 TGIs to compare, such as 30, 50, 70: the animals the study needs if every arm had each one, sized for the target power whatever the animals per arm. Scenarios, not predictions. Each adds a whole design to the run. Empty: off.'}}) tgiScenarios?: string | null,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Significance (α)', optional: true, min: '0.001', max: '0.2', initialValue: '0.05', description: 'Two-sided significance level of the comparisons with the control.'}}) alpha?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Lost before randomization', optional: true, units: '%', min: '0', max: '80', initialValue: '20', description: 'Implanted animals that never reach randomization (no take, or out of the volume window), at your facility.'}}) preLoss?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Lost after randomization', optional: true, units: '%', min: '0', max: '50', initialValue: '10', description: 'Randomized animals that cannot be analyzed (ulceration, deaths, technical failures), at your facility.'}}) postLoss?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Extra animals for losses', optional: true, choices: ['Enough for every group', 'Average loss only'], initialValue: 'Enough for every group', description: 'Average loss only needs fewer animals, but a group can end up short.'}}) lossRule?: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Confidence', optional: true, units: '%', min: '50', max: '99', initialValue: '90', description: 'How sure to be that a group keeps its animals (with Enough for every group). The protocol says what that gives for every group at once.'}}) confidence?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Volume at randomization', optional: true, units: 'mm³', min: '30', max: '500', initialValue: '100', description: 'Tumor volume when the animals are randomized and dosing starts.'}}) vRand?: number,
    // set on the study calendar, not in the form: the platform's `visible` rule, which datagrok-api's InputOptions does
    // not declare
    // @ts-expect-error
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Measurement days', optional: true, visible: 'false', initialValue: 'Mon, Wed, Fri', description: 'The weekdays the tumors are measured, e.g. Tue, Fri, also set from the study calendar.'}}) measurementDays?: string,
    // @ts-expect-error (as above)
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Randomization day', optional: true, visible: 'false', choices: ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'], initialValue: 'Monday', description: 'The day of the week the animals are randomized, measurements staying on your usual weekdays.'}}) startWeekday?: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Humane limit', optional: true, units: 'mm³', min: '500', max: '4000', initialValue: '1500', description: 'An animal is removed when its tumor passes your institution’s limit, whatever the study day. Workman 2010 advises a mean diameter of 1.5 cm at most in therapy studies, about 1,700 mm³. Your lab’s limit is set once, in the package’s settings.'}}) humane?: number,
    // set by the Plan view, which compares the schedules only while their section is open
    // @ts-expect-error (as above)
    @grok.decorators.param({type: 'bool', options: {category: 'Advanced', caption: 'Compare measurement schedules', optional: true, visible: 'false', initialValue: 'true', description: 'Also give the power by measurements per week (off: that table is empty and the plan twice as fast).'}}) compareSchedules?: boolean,
    // set by the Plan view, which shows the plan first and fills the smallest detectable effects in afterwards
    // @ts-expect-error (as above)
    @grok.decorators.param({type: 'bool', options: {category: 'Advanced', caption: 'Smallest detectable effect', optional: true, visible: 'false', initialValue: 'true', description: 'With the animals per arm set, also search each arm’s smallest detectable effect (off: left empty and the plan several times faster).'}}) detectableEffect?: boolean,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Equally good fits', optional: true, choices: ['Average their power', 'Best fit only'], initialValue: 'Average their power', description: 'For a tumor model whose published controls several fits describe equally well (Panc02, PDX-PDAC; the plan says so when they disagree): size the groups for the power averaged over those fits, or for the atlas’s best fit alone. It changes nothing for the other models.'}}) fits?: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Effects as', optional: true, choices: ['TGI', 'T/C'], description: 'How every arm’s expected effect is given: tumor growth inhibition (TGI) or the treated / control volume ratio (T/C), in percent. Empty: each arm’s own metric in the arms table.'}}) effectAs?: string,
  ): Promise<PlanFrames> {
    return planFrames(declarationOf({line, endpoint, arms, power, correction, measurementDays, n, lastDay, alpha,
      preLoss, postLoss, lossRule, confidence, vRand, startWeekday, humane, compareSchedules, detectableEffect,
      fits, effectAs}), tgisOf(tgiScenarios));
  }

  @grok.decorators.func({
    name: 'pmxAnalyzeStudy',
    description: 'Analyzes a xenograft or syngeneic efficacy study as PMX Studio’s Plan pre-specifies it, with Crux stats.tgi and stats.contrast: each treated arm against the control on the per-animal slope of log tumor volume from randomization to the evaluation day (Welch’s t test, each group its own variance, with the multiplicity correction chosen: Holm by default). The evaluation day is the first measurement day the median control reaches the endpoint volume, as in the design. Also each arm’s TGI and T/C on that day, the median time to the endpoint volume, and per animal the slope, the time to the endpoint and the mRECIST category. With Crux stats.tumor_twin, the digital twin of each treated animal: what it would measure untreated (from its volume at randomization, the study’s controls and the cell line’s model), its paired effect on the evaluation day, and each control checked against the others. Given the plan saved in Plan, the study against it: the controls against their predicted growth, the end against the predicted day, each arm’s TGI or T/C against the one assumed. Given the body weight, with Crux stats.tumor_body_weight, whether the treatments were tolerated: each group’s weight change from day 0, measured and without the tumor’s mass (its volume at the cell line’s density, or at the study’s when the excised tumors were weighed), its weight-loss grades (CTCAE v5, adapted), the animals removed early and whether each group stays within the NCI limits for a treatment (a mean loss at the nadir over 20%, or over 10% of the animals dead or removed for their weight).',
    outputs: [
      {name: 'analysis', type: 'dataframe', options: {description: 'The analysis: one row.'}},
      {name: 'comparisons', type: 'dataframe', options: {description: 'Each treated arm against the control: one row each.'}},
      {name: 'animals', type: 'dataframe', options: {description: 'Per animal: slope of log volume, time to the endpoint, mRECIST.'}},
      {name: 'volumes', type: 'dataframe', options: {description: 'Observed volumes per group and day: median and range.'}},
      {name: 'twins', type: 'dataframe', options: {description: 'Per treated animal and day: its measured volume and what it would measure untreated (median and 90% band).'}},
      {name: 'effects', type: 'dataframe', options: {description: 'Per treated animal: how much smaller it is than its untreated twin on the evaluation day, with a 90% interval.'}},
      {name: 'checks', type: 'dataframe', options: {description: 'Per control measurement: its band predicted from the other controls, and whether it falls inside.'}},
      {name: 'versusPlan', type: 'dataframe', options: {description: 'Per treated arm: the plan arm it ran as, the effect the plan assumed and the one observed with its 95% interval, and whether the assumption held.'}},
      {name: 'predicted', type: 'dataframe', options: {description: 'The controls as the plan (or, with none, the atlas) predicted them and as observed: median and range per day.'}},
      {name: 'tolerability', type: 'dataframe', options: {description: 'Per group, measured and tumor-corrected: the weight-loss grades, the mean nadir, the recoveries, the animals removed and whether it is beyond the NCI limits (empty without body weight).'}},
      {name: 'weights', type: 'dataframe', options: {description: 'Per group and day, measured and tumor-corrected: the mean body weight change from day 0 and its range (empty without body weight).'}},
    ],
  })
  static async pmxAnalyzeStudy(
    @grok.decorators.param({options: {category: 'Study', caption: 'Study', nullable: false, description: 'The measurements: one row per animal and day, day 0 the randomization (start of treatment).'}}) study: DG.DataFrame,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Animal', nullable: false, description: 'The animal (tumor) identifier, unique in the study.'}}) animal: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Group', nullable: false, description: 'The treatment group.'}}) group: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Day', nullable: false, type: 'numerical', description: 'Days from randomization.'}}) day: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Volume', nullable: false, type: 'numerical', description: 'Tumor volume (mm³).'}}) volume: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Body weight', type: 'numerical', optional: true, nullable: true, description: 'Body weight (g), weighed with the tumor (empty: no tolerability).'}}) bodyWeight: DG.Column | null,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Tumor weight', type: 'numerical', optional: true, nullable: true, description: 'The excised tumor’s weight (g), on the row of the animal’s last measurement.'}}) tumorWeight: DG.Column | null,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Removal reason', optional: true, nullable: true, description: 'Why an animal left the study (found dead, moribund, weight, tumor, ulceration, other), on any of its rows.'}}) removalReason: DG.Column | null,
    // the plan's text, not a file: the Analyze view reads the file into it
    // @ts-expect-error (the platform's `visible` rule, as above)
    @grok.decorators.param({type: 'string', options: {category: 'Study', caption: 'Plan', optional: true, nullable: true, visible: 'false', description: 'The plan saved in Plan (its .pmxplan file’s contents), to compare the study with.'}}) plan: string | null,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Cell line', nullable: false, choices: 'PmxStudio:pmxLines()', initialValue: 'MC38', description: 'The tumor model: its growth model from the atlas of historical controls predicts each treated animal untreated.'}}) line: string,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Control group', optional: true, nullable: true, description: 'The control (vehicle) group (empty: the one named vehicle, control, PBS, saline or untreated, else the first).'}}) control: string,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Endpoint volume', nullable: false, units: 'mm³', min: '500', max: '3000', initialValue: '1000', description: 'The effects are read on the first measurement day the median control reaches this volume, as the design defines the end of the study.'}}) endpoint: number,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Multiple comparisons', nullable: false, choices: ['holm', 'bonferroni', 'none'], initialValue: 'holm', description: 'Correction over the arms compared with the control, each by Welch’s t test: Holm (recommended), Bonferroni, or none (each comparison at alpha, as in an exploratory study).'}}) correction: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Evaluation day', optional: true, nullable: true, description: 'Leave empty for the day the median control reaches the endpoint volume (the last day if it never does).'}}) evaluationDay?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Significance (α)', optional: true, min: '0.001', max: '0.2', initialValue: '0.05', description: 'Two-sided significance level of the comparisons with the control.'}}) alpha?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Humane limit, weight loss', optional: true, units: '%', min: '5', max: '50', initialValue: '20', description: 'An animal is removed when it loses this much weight (your institution’s limit).'}}) weightLimit?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Humane limit, tumor', optional: true, units: 'mm³', min: '500', max: '4000', initialValue: '1500', description: 'An animal is removed when its tumor passes your institution’s limit. Workman 2010 advises a mean diameter of 1.5 cm at most in therapy studies, about 1,700 mm³. Your lab’s limit is set once, in the package’s settings.'}}) humane?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Tumor density', optional: true, nullable: true, units: 'mg/mm³', min: '0.1', max: '3', description: 'The tumor’s mass per mm³ of caliper volume (empty: the cell line’s typical one).'}}) density?: number,
    // set by the Analyze view, which shows the analysis first and fills the twins in afterwards
    // @ts-expect-error (the platform's `visible` rule, as above)
    @grok.decorators.param({type: 'bool', options: {category: 'Advanced', caption: 'Digital twin', optional: true, visible: 'false', initialValue: 'true', description: 'Also compare each treated animal with its untreated twin (off: those tables are empty and the analysis several times faster).'}}) withTwin?: boolean,
  ): Promise<AnalysisFrames> {
    return analysisFrames(analysisDeclarationOf({study, animal, group, day, volume, bodyWeight, tumorWeight,
      removalReason, plan, line, control, endpoint, correction, evaluationDay, alpha, weightLimit, humane, density,
      withTwin}));
  }
}
