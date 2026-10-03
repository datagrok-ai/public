/** Cross-model run comparison (F-29): a criterion's candidate runs, and Compute2's Compare Runs.
 * The plan is pure (src/model/compare-runs.ts); this reads the runs and opens the tool. */
import * as DG from 'datagrok-api/dg';
import {CriterionRuns} from './model-runs';
import {COMPARE_RUNS_FUNCTION, CompareCandidate, comparePlan, ComparePlan} from './model/compare-runs';

/** Compute2's Compare Runs, or null when the package is not installed (or not shared with the user). */
export function compareRunsTool(): DG.Func | null {
  const [pkg, name] = COMPARE_RUNS_FUNCTION.split(':');
  return DG.Func.find({package: pkg, name})[0] ?? null;
}

/** The criterion's candidates as the plan lists them, newest first. A run counts as another
 * user's when its recorded author is not the current user. */
export function compareCandidates(runs: CriterionRuns, me: string = DG.User.current().id): CompareCandidate[] {
  return runs.candidates.map((i) => ({id: i.run.id, runId: i.run.run_id ?? null,
    title: i.run.title ?? i.run.run_id ?? i.run.model_name, model: i.run.model_name,
    modelName: i.model?.name ?? i.run.model_name, on: i.run.created_on ? i.run.created_on.toISOString() : null,
    selected: i.run.id === runs.selectedRunId, byOther: !!i.run.run_author && i.run.run_author !== me}));
}

export function criterionComparePlan(runs: CriterionRuns): ComparePlan {
  return comparePlan(compareCandidates(runs), compareRunsTool() !== null);
}

/** Opens Compare Runs in a view of its own, empty: it takes no runs. */
export async function openCompareRuns(): Promise<void> {
  const f = compareRunsTool();
  if (!f)
    throw new Error(`${COMPARE_RUNS_FUNCTION} is not installed on this server`);
  await f.apply({});
}

/** The plan as a table: one row per candidate. */
export function comparePlanFrame(p: ComparePlan, name: string): DG.DataFrame {
  const runs = p.groups.flatMap((g) => g.runs);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Model', runs.map((r) => r.modelName)),
    DG.Column.fromStrings('Run', runs.map((r) => r.title)),
    DG.Column.fromStrings('Run id', runs.map((r) => r.runId ?? '')),
    DG.Column.fromStrings('Saved on', runs.map((r) => r.on?.slice(0, 10) ?? '')),
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'Selected', runs.map((r) => r.selected)),
  ]);
  df.name = name;
  return df;
}
