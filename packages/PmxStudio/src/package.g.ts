import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: PMX Studio
//description: Plans a xenograft or syngeneic efficacy study de novo from the atlas of historical controls.
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.role: app
export async function pmxStudioApp(path?: string) : Promise<any> {
  return await PackageFunctions.pmxStudioApp(path);
}

//description: Opens a PMX Studio plan in Plan.
//input: string content 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: pmxplan
export async function pmxPlanFileHandler(content: string) : Promise<any> {
  return await PackageFunctions.pmxPlanFileHandler(content);
}

//description: Control-growth priors of the PMX Studio atlas (reduced Gompertz) per syngeneic cell line: study-level initial growth rate alpha (median and 95% interval of a new study), deceleration beta, between-study tau, effective sample size, non-take rate, papers and animals behind each prior.
//output: dataframe result
export async function pmxPriors() : Promise<any> {
  return await PackageFunctions.pmxPriors();
}

//description: Recomputes the MAP prior of a new study initial growth rate (log alpha) for a cell line with Crux stats.map_prior from the atlas per-paper estimates: one row per mixture component.
//input: string line { choices: ["MC38","CT26","B16F10"]; description: Syngeneic cell line. }
//output: dataframe result
export async function pmxMapPrior(line: string) : Promise<any> {
  return await PackageFunctions.pmxMapPrior(line);
}

//description: Designs a xenograft or syngeneic efficacy study de novo with Crux stats.tumor_study_design: from the cell line, the endpoint volume and the treated arms (each with its kind of therapy and expected TGI or T/C), predicts the end of the study (when the median control reaches the endpoint) and its uncertainty, each arm’s kill rate, the analysable animals per arm for the declared assurance (power over the historical-control prior) with the multiplicity correction chosen (Dunnett by default), the animals to implant for them, the power curves and the predicted control and treated arms. Given a number of animals per arm, the power there and each arm’s minimum detectable effect. Given an end of study, the study is evaluated that day; the power by end of study shows what ending earlier or later costs, and the power by measurement schedule what measuring less often costs.
//input: string line = 'MC38' { category: Tumor model; caption: Cell line; choices: ["MC38","CT26","B16F10"]; description: Syngeneic cell line: its growth prior comes from the atlas of historical controls. }
//input: double endpoint = 1500 { category: Tumor model; caption: Endpoint volume; units: mm³; min: 500; max: 3000; description: The study ends when the median control reaches it. }
//input: dataframe arms { category: Treatment arms; caption: Arms; optional: true; description: The treated arms, each compared with the control: name, therapy (small_molecule, biologic, immunotherapy, immune_agonist, cell_therapy, vaccine_peptide, oligo_rna), metric (TGI or TC) and the expected effect in percent. Empty: one small-molecule arm at TGI 60%. }
//input: double power = 80 { category: Design; caption: Target power; units: %; min: 50; max: 99; description: Chance of detecting each arm’s effect. }
//input: string correction = 'dunnett' { category: Design; caption: Multiple comparisons; choices: ["dunnett","bonferroni","none"]; description: Correction over the arms compared with the control: Dunnett (recommended), Bonferroni, or none (each comparison at alpha, as in an exploratory study). }
//input: int n { category: Design; caption: Animals per arm; optional: true; nullable: true; min: 2; max: 100; description: Empty for the number the target power needs, or set it to see the power and smallest detectable effect with that many. }
//input: int lastDay { category: Design; caption: End of study; optional: true; nullable: true; units: day; min: 5; max: 120; description: Leave empty to end when the median control reaches the endpoint. }
//input: double alpha = 0.05 { category: Advanced; caption: Significance (α); optional: true; min: 0.001; max: 0.2; description: Two-sided significance level of the comparisons with the control. }
//input: double preLoss = 20 { category: Advanced; caption: Lost before randomization; optional: true; units: %; min: 0; max: 80; description: Implanted animals that never reach randomization (no take, or out of the volume window), at your facility. }
//input: double postLoss = 10 { category: Advanced; caption: Lost after randomization; optional: true; units: %; min: 0; max: 50; description: Randomized animals that cannot be analyzed (ulceration, deaths, technical failures), at your facility. }
//input: string lossRule = 'Enough for every group' { category: Advanced; caption: Extra animals for losses; optional: true; choices: ["Enough for every group","Average loss only"]; description: Average loss only needs fewer animals, but a group can end up short. }
//input: double confidence = 90 { category: Advanced; caption: Confidence; optional: true; units: %; min: 50; max: 99; description: How sure to be that every group keeps its animals (with Enough for every group). }
//input: double vRand = 100 { category: Advanced; caption: Volume at randomization; optional: true; units: mm³; min: 30; max: 500; description: Tumor volume when the animals are randomized and dosing starts. }
//input: string measurementDays = 'Mon, Wed, Fri' { category: Advanced; caption: Measurement days; optional: true; visible: false; description: The weekdays the tumors are measured, e.g. Tue, Fri, also set from the study calendar. }
//input: string startWeekday = 'Monday' { category: Advanced; caption: Randomization day; optional: true; visible: false; choices: ["Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"]; description: The day of the week the animals are randomized, measurements staying on your usual weekdays. }
//input: double humane = 2000 { category: Advanced; caption: Humane limit; optional: true; units: mm³; min: 500; max: 4000; description: An animal is removed when its tumor passes your institution’s limit, whatever the study day. }
//input: bool compareSchedules = true { category: Advanced; caption: Compare measurement schedules; optional: true; visible: false; description: Also give the power by measurements per week (off: that table is empty and the plan twice as fast). }
//output: dataframe plan { description: The study: one row. }
//output: dataframe armResults { description: The treated arms: one row each. }
//output: dataframe powerCurve { description: Assurance (power averaged over the labs of the atlas) per arm and number of analysable animals per arm. }
//output: dataframe curves { description: Predicted volumes of the control and each treated arm, day by day from randomization: median and 50% / 90% bands. }
//output: dataframe duration { description: Each arm’s power by the day the study ends, with the share of controls past the endpoint and over the humane limit. }
//output: dataframe schedules { description: Each arm’s power and the analyzable animals the study needs, by measurements per week (empty when not asked for). }
//output: dataframe calendar { description: The study day by day from randomization: what is done each day. }
export async function pmxPlanStudy(line: string, endpoint: number, arms: DG.DataFrame, power: number, correction: string, n?: number, lastDay?: number, alpha?: number, preLoss?: number, postLoss?: number, lossRule?: string, confidence?: number, vRand?: number, measurementDays?: string, startWeekday?: string, humane?: number, compareSchedules?: boolean) : Promise<any> {
  return await PackageFunctions.pmxPlanStudy(line, endpoint, arms, power, correction, n, lastDay, alpha, preLoss, postLoss, lossRule, confidence, vRand, measurementDays, startWeekday, humane, compareSchedules);
}

//description: Analyzes a xenograft or syngeneic efficacy study as PMX Studio’s Plan pre-specifies it, with Crux stats.tgi and stats.contrast: each treated arm against the control on the per-animal slope of log tumor volume from randomization to the evaluation day (two-sample t test, with the multiplicity correction chosen: Dunnett by default). The evaluation day is the first measurement day the median control reaches the endpoint volume, as in the design. Also each arm’s TGI and T/C on that day, the median time to the endpoint volume, and per animal the slope, the time to the endpoint and the mRECIST category. With Crux stats.tumor_twin, the digital twin of each treated animal: what it would measure untreated (from its volume at randomization, the study’s controls and the cell line’s model), its paired effect on the evaluation day, and each control checked against the others. Given the plan saved in Plan, the study against it: the controls against their predicted growth, the end against the predicted day, each arm’s TGI or T/C against the one assumed. Given the body weight, with Crux stats.tumor_body_weight, whether the treatments were tolerated: each group’s weight change from day 0, measured and without the tumor’s mass (its volume at the cell line’s density, or at the study’s when the excised tumors were weighed), its weight-loss grades (CTCAE v5, adapted), the animals removed early and whether each group stays within the NCI limits for a treatment (a mean loss at the nadir over 20%, or over 10% of the animals dead or removed for their weight).
//input: dataframe study { category: Study; caption: Study; nullable: false; description: The measurements: one row per animal and day, day 0 the randomization (start of treatment). }
//input: column animal { category: Study; caption: Animal; description: The animal (tumor) identifier, unique in the study. }
//input: column group { category: Study; caption: Group; description: The treatment group. }
//input: column day { category: Study; caption: Day; type: numerical; description: Days from randomization. }
//input: column volume { category: Study; caption: Volume; type: numerical; description: Tumor volume (mm³). }
//input: column bodyWeight { category: Study; caption: Body weight; type: numerical; optional: true; nullable: true; description: Body weight (g), weighed with the tumor (empty: no tolerability). }
//input: column tumorWeight { category: Study; caption: Tumor weight; type: numerical; optional: true; nullable: true; description: The excised tumor’s weight (g), on the row of the animal’s last measurement. }
//input: column removalReason { category: Study; caption: Removal reason; optional: true; nullable: true; description: Why an animal left the study (found dead, moribund, weight, tumor, ulceration, other), on any of its rows. }
//input: string plan { category: Study; caption: Plan; optional: true; nullable: true; visible: false; description: The plan saved in Plan (its .pmxplan file’s contents), to compare the study with. }
//input: string line = 'MC38' { category: Analysis; caption: Cell line; choices: ["MC38","CT26","B16F10"]; description: The tumor model: its growth model from the atlas of historical controls predicts each treated animal untreated. }
//input: string control { category: Analysis; caption: Control group; optional: true; nullable: true; description: The control (vehicle) group (empty: the one named vehicle, control, PBS, saline or untreated, else the first). }
//input: double endpoint = 1500 { category: Analysis; caption: Endpoint volume; units: mm³; min: 500; max: 3000; description: The effects are read on the first measurement day the median control reaches this volume, as the design defines the end of the study. }
//input: string correction = 'dunnett' { category: Analysis; caption: Multiple comparisons; choices: ["dunnett","bonferroni","none"]; description: Correction over the arms compared with the control: Dunnett (recommended), Bonferroni, or none (each comparison at alpha, as in an exploratory study). }
//input: double evaluationDay { category: Advanced; caption: Evaluation day; optional: true; nullable: true; description: Leave empty for the day the median control reaches the endpoint volume (the last day if it never does). }
//input: double alpha = 0.05 { category: Advanced; caption: Significance (α); optional: true; min: 0.001; max: 0.2; description: Two-sided significance level of the comparisons with the control. }
//input: double weightLimit = 20 { category: Advanced; caption: Humane limit, weight loss; optional: true; units: %; min: 5; max: 50; description: An animal is removed when it loses this much weight (your institution’s limit). }
//input: double humane = 2000 { category: Advanced; caption: Humane limit, tumor; optional: true; units: mm³; min: 500; max: 4000; description: An animal is removed when its tumor passes your institution’s limit. }
//input: double density { category: Advanced; caption: Tumor density; optional: true; nullable: true; units: mg/mm³; min: 0.1; max: 3; description: The tumor’s mass per mm³ of caliper volume (empty: the cell line’s typical one). }
//output: dataframe analysis { description: The analysis: one row. }
//output: dataframe comparisons { description: Each treated arm against the control: one row each. }
//output: dataframe animals { description: Per animal: slope of log volume, time to the endpoint, mRECIST. }
//output: dataframe volumes { description: Observed volumes per group and day: median and range. }
//output: dataframe twins { description: Per treated animal and day: its measured volume and what it would measure untreated (median and 90% band). }
//output: dataframe effects { description: Per treated animal: how much smaller it is than its untreated twin on the evaluation day, with a 90% interval. }
//output: dataframe checks { description: Per control measurement: its band predicted from the other controls, and whether it falls inside. }
//output: dataframe versusPlan { description: Per treated arm: the plan arm it ran as, the effect the plan assumed and the one observed with its 95% interval, and whether the assumption held. }
//output: dataframe predicted { description: The controls as the plan (or, with none, the atlas) predicted them and as observed: median and range per day. }
//output: dataframe tolerability { description: Per group, measured and tumor-corrected: the weight-loss grades, the mean nadir, the recoveries, the animals removed and whether it is beyond the NCI limits (empty without body weight). }
//output: dataframe weights { description: Per group and day, measured and tumor-corrected: the mean body weight change from day 0 and its range (empty without body weight). }
export async function pmxAnalyzeStudy(study: DG.DataFrame, animal: DG.Column, group: DG.Column, day: DG.Column, volume: DG.Column, bodyWeight: any, tumorWeight: any, removalReason: any, plan: any, line: string, control: string, endpoint: number, correction: string, evaluationDay?: number, alpha?: number, weightLimit?: number, humane?: number, density?: number) : Promise<any> {
  return await PackageFunctions.pmxAnalyzeStudy(study, animal, group, day, volume, bodyWeight, tumorWeight, removalReason, plan, line, control, endpoint, correction, evaluationDay, alpha, weightLimit, humane, density);
}
