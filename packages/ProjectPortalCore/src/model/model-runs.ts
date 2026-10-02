/** A model run saved to a program: what it is annotated with, and where in the workspace it goes.
 * Pure data, so the rules hold wherever a run is read back. */

/** The FuncCall options a portal run carries. Options persist with the saved call and can be
 * queried (`options.portalProgram = "EXB-203"`), so a run is findable from the call alone,
 * even when no `model_run` row could be written. */
export const RUN_OPTIONS = {
  portfolio: 'portalPortfolio',
  program: 'portalProgram',
  stage: 'portalStage',
  unitOperation: 'portalUnitOperation',
  runType: 'portalRunType',
  purpose: 'portalPurpose',
  /** Ids of the tables uploaded for the call's dataframes. */
  tables: 'portalTables',
} as const;

export const RUN_PURPOSES = ['screening', 'characterization', 'scale-up', 'troubleshooting', 'verification'] as const;
export type RunPurpose = typeof RUN_PURPOSES[number];

export interface RunAnnotations {
  portfolio: string;
  program: string;
  stage?: string | null;
  unitOperation?: string | null;
  runType?: string | null;
  purpose?: RunPurpose | null;
  title: string;
  description?: string | null;
}

/** The options to write on the call: the portal keys, and the `title`, `description` and `tags`
 * the run history of Compute shows. Tags already on the call are kept. */
export function runOptions(a: RunAnnotations, existingTags: string[] = []): Record<string, any> {
  const options: Record<string, any> = {
    [RUN_OPTIONS.portfolio]: a.portfolio,
    [RUN_OPTIONS.program]: a.program,
    title: a.title,
  };
  const optional: [string, string | null | undefined][] = [
    [RUN_OPTIONS.stage, a.stage], [RUN_OPTIONS.unitOperation, a.unitOperation],
    [RUN_OPTIONS.runType, a.runType], [RUN_OPTIONS.purpose, a.purpose], ['description', a.description]];
  for (const [k, v] of optional) {
    if (v)
      options[k] = v;
  }
  const tags = [...existingTags];
  for (const t of ['portal', a.program, a.stage])
    if (t && !tags.includes(t))
      tags.push(t);
  options.tags = tags;
  return options;
}

/** Reads the portal annotations back from a call's options; null when the call is not a portal run. */
export function readRunAnnotations(options: Record<string, any> | null | undefined): RunAnnotations | null {
  if (!options || !options[RUN_OPTIONS.program])
    return null;
  const s = (k: string) => options[k] == null ? null : String(options[k]);
  return {
    portfolio: s(RUN_OPTIONS.portfolio) ?? '',
    program: s(RUN_OPTIONS.program)!,
    stage: s(RUN_OPTIONS.stage),
    unitOperation: s(RUN_OPTIONS.unitOperation),
    runType: s(RUN_OPTIONS.runType),
    purpose: s(RUN_OPTIONS.purpose) as RunPurpose | null,
    title: s('title') ?? '',
    description: s('description'),
  };
}

/** Where under the program's folder a run is linked: the folder asked for, else the first folder
 * the program's current stage declares (so a run saved at process lock lands in
 * Drug Substance/Process Lock), else the program folder itself ([]). */
export function runFolderPath(stageFolders: string | null | undefined, asked?: string | null): string[] {
  const split = (p: string) => p.split('/').map((s) => s.trim()).filter((s) => s.length > 0);
  if (asked && split(asked).length > 0)
    return split(asked);
  const first = (stageFolders ?? '').split('|').map((f) => f.trim()).find((f) => f.length > 0);
  return first ? split(first) : [];
}
