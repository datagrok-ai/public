/** Cross-model run comparison (F-29, FR-55): a criterion's candidate runs, handed to Compute2's
 * Compare Runs tool. Not a portal build: the tool does the matching of inputs and outputs across
 * models. Compare Runs takes no runs as input today, so the portal cannot seed it; it opens the
 * tool and lists what to add, model by model. Pure: the caller reads the runs and says whether the
 * tool is installed. */

/** Compute2's run comparison tool (Model Hub → Tools → Compare Runs). */
export const COMPARE_RUNS_FUNCTION = 'Compute2:CompareRuns';

/** A candidate run as the plan lists it. */
export interface CompareCandidate {
  /** The model_run id. */
  id: string;
  /** The saved call's id; null when the run has none (Compare Runs cannot add it). */
  runId: string | null;
  title: string;
  /** The model's function, as nqName. */
  model: string;
  /** The model's name in the catalogue, else its function. */
  modelName: string;
  /** ISO time. */
  on: string | null;
  selected: boolean;
}

export interface CompareGroup {
  model: string;
  modelName: string;
  runs: CompareCandidate[];
}

export interface ComparePlan {
  /** One group per model, in the order the models first appear (newest run first). */
  groups: CompareGroup[];
  runs: number;
  models: number;
  /** Runs of more than one model: Compare Runs matches them by parameter name and units. */
  crossModel: boolean;
  /** Why Compare Runs cannot be opened for them; null when it can. */
  unavailable: string | null;
  /** What to do in the tool, one step a line. */
  steps: string[];
  /** Runs without a saved call: listed, but the tool cannot add them. */
  notAddable: number;
}

/** The candidates grouped by model, and what to do in Compare Runs. [toolInstalled] is whether
 * COMPARE_RUNS_FUNCTION exists on this server. */
export function comparePlan(candidates: CompareCandidate[], toolInstalled: boolean): ComparePlan {
  const groups: CompareGroup[] = [];
  for (const c of candidates) {
    let g = groups.find((x) => x.model === c.model);
    if (!g)
      groups.push(g = {model: c.model, modelName: c.modelName, runs: []});
    g.runs.push(c);
  }
  const addable = candidates.filter((c) => c.runId);
  const unavailable = !toolInstalled ? 'Compute2\'s Compare Runs is not installed on this server' :
    addable.length < 2 ? `Comparing needs two runs with a saved call; ${addable.length === 1 ? 'one' : 'none'} here` : null;
  const steps = unavailable ? [] : [
    'Open Compare Runs',
    ...groups.filter((g) => g.runs.some((r) => r.runId)).map((g) =>
      `Pick ${g.modelName}, tick ${g.runs.filter((r) => r.runId).map((r) => `"${r.title}"`).join(', ')}, then Add selected runs`),
    ...(groups.length > 1 ? ['Runs of different models are matched by input and output name and units; check the pairings it proposes'] : []),
  ];
  return {groups, runs: candidates.length, models: groups.length, crossModel: groups.length > 1, unavailable, steps,
    notAddable: candidates.length - addable.length};
}

/** "3 runs of 2 models" */
export function describePlan(p: ComparePlan): string {
  const n = (k: number, one: string) => `${k} ${one}${k === 1 ? '' : 's'}`;
  return p.runs === 0 ? 'No candidate runs' : `${n(p.runs, 'run')} of ${n(p.models, 'model')}`;
}
