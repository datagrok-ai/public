/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PlanFrames, declarationOf, mapPriorFrame, planFrames, priorsFrame} from './frames';
import {openPlanFile, planView} from './plan-view';
import {AnalysisFrames, analysisDeclarationOf, analysisFrames} from './analyze';
import {analyzeView} from './analyze-view';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.app({name: 'PMX Studio',
    description: 'Plans a xenograft or syngeneic efficacy study de novo from the atlas of historical controls.'})
  static async pmxStudioApp(
    // the app's URL path, which the platform passes: /analyze opens Analyze, anything else Plan
    @grok.decorators.param({options: {metaUrl: true, optional: true}}) path?: string): Promise<DG.ViewBase> {
    // ponytail: the stage from the URL until the four-stage hub exists
    return path?.replace(/^\//, '').startsWith('analyze') ? analyzeView() : planView();
  }

  // the extension as a literal (the decorators' generator reads no constants): plan-file.ts's PLAN_EXT
  @grok.decorators.fileHandler({ext: 'pmxplan', description: 'Opens a PMX Studio plan in Plan.'})
  static async pmxPlanFileHandler(content: string): Promise<DG.DataFrame[]> {
    return openPlanFile(content);
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
    @grok.decorators.param({options: {choices: ['MC38', 'CT26', 'B16F10'], description: 'Syngeneic cell line.'}}) line: string): Promise<DG.DataFrame> {
    return mapPriorFrame(line);
  }

  @grok.decorators.func({
    name: 'pmxPlanStudy',
    description: 'Designs a xenograft or syngeneic efficacy study de novo with Crux stats.tumor_study_design: from the cell line, the endpoint volume and the treated arms (each with its kind of therapy and expected TGI or T/C), predicts the end of the study (when the median control reaches the endpoint) and its uncertainty, each arm’s kill rate, the analysable animals per arm for the declared assurance (power over the historical-control prior) with the multiplicity correction chosen (Dunnett by default), the animals to implant for them, the power curves and the predicted control and treated arms. Given a number of animals per arm, the power there and each arm’s minimum detectable effect. Given an end of study, the study is evaluated that day; the power by end of study shows what ending earlier or later costs, and the power by measurement schedule what measuring less often costs.',
    outputs: [
      {name: 'plan', type: 'dataframe', options: {description: 'The study: one row.'}},
      {name: 'armResults', type: 'dataframe', options: {description: 'The treated arms: one row each.'}},
      {name: 'powerCurve', type: 'dataframe', options: {description: 'Assurance (power averaged over the labs of the atlas) per arm and number of analysable animals per arm.'}},
      {name: 'curves', type: 'dataframe', options: {description: 'Predicted volumes of the control and each treated arm, day by day from randomization: median and 50% / 90% bands.'}},
      {name: 'duration', type: 'dataframe', options: {description: 'Each arm’s power by the day the study ends, with the share of controls past the endpoint and over the humane limit.'}},
      {name: 'schedules', type: 'dataframe', options: {description: 'Each arm’s power and the analyzable animals the study needs, by measurements per week (empty when not asked for).'}},
      {name: 'calendar', type: 'dataframe', options: {description: 'The study day by day from randomization: what is done each day.'}},
    ],
  })
  static async pmxPlanStudy(
    @grok.decorators.param({options: {category: 'Tumor model', caption: 'Cell line', choices: ['MC38', 'CT26', 'B16F10'], initialValue: 'MC38', description: 'Syngeneic cell line: its growth prior comes from the atlas of historical controls.'}}) line: string,
    @grok.decorators.param({options: {category: 'Tumor model', caption: 'Endpoint volume', units: 'mm³', min: '500', max: '3000', initialValue: '1500', description: 'The study ends when the median control reaches it.'}}) endpoint: number,
    @grok.decorators.param({options: {category: 'Treatment arms', caption: 'Arms', optional: true, description: 'The treated arms, each compared with the control: name, therapy (small_molecule, biologic, immunotherapy, immune_agonist, cell_therapy, vaccine_peptide, oligo_rna), metric (TGI or TC) and the expected effect in percent. Empty: one small-molecule arm at TGI 60%.'}}) arms: DG.DataFrame,
    @grok.decorators.param({options: {category: 'Design', caption: 'Target power', units: '%', min: '50', max: '99', initialValue: '80', description: 'Chance of detecting each arm’s effect.'}}) power: number,
    @grok.decorators.param({options: {category: 'Design', caption: 'Multiple comparisons', choices: ['dunnett', 'bonferroni', 'none'], initialValue: 'dunnett', description: 'Correction over the arms compared with the control: Dunnett (recommended), Bonferroni, or none (each comparison at alpha, as in an exploratory study).'}}) correction: string,
    @grok.decorators.param({type: 'int', options: {category: 'Design', caption: 'Animals per arm', optional: true, nullable: true, min: '2', max: '100', description: 'Empty for the number the target power needs, or set it to see the power and smallest detectable effect with that many.'}}) n?: number,
    @grok.decorators.param({type: 'int', options: {category: 'Design', caption: 'End of study', optional: true, nullable: true, units: 'day', min: '5', max: '120', description: 'Leave empty to end when the median control reaches the endpoint.'}}) lastDay?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Significance (α)', optional: true, min: '0.001', max: '0.2', initialValue: '0.05', description: 'Two-sided significance level of the comparisons with the control.'}}) alpha?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Lost before randomization', optional: true, units: '%', min: '0', max: '80', initialValue: '20', description: 'Implanted animals that never reach randomization (no take, or out of the volume window), at your facility.'}}) preLoss?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Lost after randomization', optional: true, units: '%', min: '0', max: '50', initialValue: '10', description: 'Randomized animals that cannot be analyzed (ulceration, deaths, technical failures), at your facility.'}}) postLoss?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Extra animals for losses', optional: true, choices: ['Enough for every group', 'Average loss only'], initialValue: 'Enough for every group', description: 'Average loss only needs fewer animals, but a group can end up short.'}}) lossRule?: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Confidence', optional: true, units: '%', min: '50', max: '99', initialValue: '90', description: 'How sure to be that every group keeps its animals (with Enough for every group).'}}) confidence?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Volume at randomization', optional: true, units: 'mm³', min: '30', max: '500', initialValue: '100', description: 'Tumor volume when the animals are randomized and dosing starts.'}}) vRand?: number,
    // set on the study calendar, not in the form: the platform's `visible` rule, which datagrok-api's InputOptions does
    // not declare
    // @ts-expect-error
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Measurement days', optional: true, visible: 'false', initialValue: 'Mon, Wed, Fri', description: 'The weekdays the tumors are measured, e.g. Tue, Fri, also set from the study calendar.'}}) measurementDays?: string,
    // @ts-expect-error (as above)
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Randomization day', optional: true, visible: 'false', choices: ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'], initialValue: 'Monday', description: 'The day of the week the animals are randomized, measurements staying on your usual weekdays.'}}) startWeekday?: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Humane limit', optional: true, units: 'mm³', min: '500', max: '4000', initialValue: '2000', description: 'An animal is removed when its tumor passes your institution’s limit, whatever the study day.'}}) humane?: number,
    // set by the Plan view, which compares the schedules only while their section is open
    // @ts-expect-error (as above)
    @grok.decorators.param({type: 'bool', options: {category: 'Advanced', caption: 'Compare measurement schedules', optional: true, visible: 'false', initialValue: 'true', description: 'Also give the power by measurements per week (off: that table is empty and the plan twice as fast).'}}) compareSchedules?: boolean,
  ): Promise<PlanFrames> {
    return planFrames(declarationOf({line, endpoint, arms, power, correction, measurementDays, n, lastDay, alpha,
      preLoss, postLoss, lossRule, confidence, vRand, startWeekday, humane, compareSchedules}));
  }

  @grok.decorators.func({
    name: 'pmxAnalyzeStudy',
    description: 'Analyzes a xenograft or syngeneic efficacy study as PMX Studio’s Plan pre-specifies it, with Crux stats.tgi and stats.contrast: each treated arm against the control on the per-animal slope of log tumor volume from randomization to the evaluation day (two-sample t test, with the multiplicity correction chosen: Dunnett by default). The evaluation day is the first measurement day the median control reaches the endpoint volume, as in the design. Also each arm’s TGI and T/C on that day, the median time to the endpoint volume, and per animal the slope, the time to the endpoint and the mRECIST category. With Crux stats.tumor_twin, the digital twin of each treated animal: what it would measure untreated (from its volume at randomization, the study’s controls and the cell line’s model), its paired effect on the evaluation day, and each control checked against the others. Given the plan saved in Plan, the study against it: the controls against their predicted growth, the end against the predicted day, each arm’s TGI or T/C against the one assumed. Given the body weight, with Crux stats.tumor_body_weight, whether the treatments were tolerated: each group’s weight change from day 0, measured and without the tumor’s mass (its volume at the cell line’s density, or at the study’s when the excised tumors were weighed), its weight-loss grades (CTCAE v5, adapted), the animals removed early and whether each group stays within the NCI limits for a treatment (a mean loss at the nadir over 20%, or over 10% of the animals dead or removed for their weight).',
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
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Animal', description: 'The animal (tumor) identifier, unique in the study.'}}) animal: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Group', description: 'The treatment group.'}}) group: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Day', type: 'numerical', description: 'Days from randomization.'}}) day: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Volume', type: 'numerical', description: 'Tumor volume (mm³).'}}) volume: DG.Column,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Body weight', type: 'numerical', optional: true, nullable: true, description: 'Body weight (g), weighed with the tumor (empty: no tolerability).'}}) bodyWeight: DG.Column | null,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Tumor weight', type: 'numerical', optional: true, nullable: true, description: 'The excised tumor’s weight (g), on the row of the animal’s last measurement.'}}) tumorWeight: DG.Column | null,
    @grok.decorators.param({type: 'column', options: {category: 'Study', caption: 'Removal reason', optional: true, nullable: true, description: 'Why an animal left the study (found dead, moribund, weight, tumor, ulceration, other), on any of its rows.'}}) removalReason: DG.Column | null,
    // the plan's text, not a file: a saved run then carries its plan (u2 funcForm saves every run, and a file
    // param keeps it from saving); the Analyze view reads the file into it
    // @ts-expect-error (the platform's `visible` rule, as above)
    @grok.decorators.param({type: 'string', options: {category: 'Study', caption: 'Plan', optional: true, nullable: true, visible: 'false', description: 'The plan saved in Plan (its .pmxplan file’s contents), to compare the study with.'}}) plan: string | null,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Cell line', choices: ['MC38', 'CT26', 'B16F10'], initialValue: 'MC38', description: 'The tumor model: its growth model from the atlas of historical controls predicts each treated animal untreated.'}}) line: string,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Control group', optional: true, nullable: true, description: 'The control (vehicle) group (empty: the one named vehicle, control, PBS, saline or untreated, else the first).'}}) control: string,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Endpoint volume', units: 'mm³', min: '500', max: '3000', initialValue: '1500', description: 'The effects are read on the first measurement day the median control reaches this volume, as the design defines the end of the study.'}}) endpoint: number,
    @grok.decorators.param({options: {category: 'Analysis', caption: 'Multiple comparisons', choices: ['dunnett', 'bonferroni', 'none'], initialValue: 'dunnett', description: 'Correction over the arms compared with the control: Dunnett (recommended), Bonferroni, or none (each comparison at alpha, as in an exploratory study).'}}) correction: string,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Evaluation day', optional: true, nullable: true, description: 'Leave empty for the day the median control reaches the endpoint volume (the last day if it never does).'}}) evaluationDay?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Significance (α)', optional: true, min: '0.001', max: '0.2', initialValue: '0.05', description: 'Two-sided significance level of the comparisons with the control.'}}) alpha?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Humane limit, weight loss', optional: true, units: '%', min: '5', max: '50', initialValue: '20', description: 'An animal is removed when it loses this much weight (your institution’s limit).'}}) weightLimit?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Humane limit, tumor', optional: true, units: 'mm³', min: '500', max: '4000', initialValue: '2000', description: 'An animal is removed when its tumor passes your institution’s limit.'}}) humane?: number,
    @grok.decorators.param({options: {category: 'Advanced', caption: 'Tumor density', optional: true, nullable: true, units: 'mg/mm³', min: '0.1', max: '3', description: 'The tumor’s mass per mm³ of caliper volume (empty: the cell line’s typical one).'}}) density?: number,
  ): Promise<AnalysisFrames> {
    return analysisFrames(analysisDeclarationOf({study, animal, group, day, volume, bodyWeight, tumorWeight,
      removalReason, plan, line, control, endpoint, correction, evaluationDay, alpha, weightLimit, humane, density}));
  }
}
