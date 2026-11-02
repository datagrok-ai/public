import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {SaveRunDialog, SaveRunResult} from './SaveRunDialog';
import {Place, addRunsToSpace, addToSpace, placeName} from './run-storage';
import {showVueDialog} from '../../utils';

/** Links several runs into a space under their own names; runs already there are skipped. */
export async function addSeveralToSpace(func: DG.Func, runs: DG.FuncCall[]): Promise<Place | null> {
  const result = await showVueDialog<SaveRunResult>('Add to space', SaveRunDialog,
    {func, linkedRuns: runs.map((run) => Vue.markRaw(run))}, {width: 760, height: 520});
  if (!result?.place || !result.outcomes)
    return null;
  const where = placeName(result.place);
  const {added, failed} = await addRunsToSpace(result.place, result.outcomes);
  const skipped = result.outcomes.filter((outcome) => outcome.kind === 'skip').length;
  if (added.length > 0) {
    grok.shell.info(`Added ${added.length} ${added.length > 1 ? 'runs' : 'run'} to ${where}` +
      (skipped ? ` (${skipped} already there)` : ''));
  }
  if (failed.length > 0) {
    grok.shell.error(`Couldn't add ${failed.map(({run}) => getRunName(run)).join(', ')} to ${where}: ` +
      failed[0].error);
  }
  return added.length > 0 ? result.place : null;
}

/** Links an existing run into a space chosen in the Save dialog's tree, replacing a run with the same name there. */
export async function addRunToSpace(func: DG.Func, run: DG.FuncCall): Promise<Place | null> {
  const result = await showVueDialog<SaveRunResult>('Add to space', SaveRunDialog,
    {func, name: getRunName(run), linkedRun: Vue.markRaw(run)}, {width: 760, height: 480});
  if (!result?.place)
    return null;
  try {
    await addToSpace(run, result.place, result.metadata.name, result.target.replace);
    grok.shell.info(`Added to ${placeName(result.place)}`);
    return result.place;
  } catch (e: any) {
    grok.shell.error(`Couldn't add the run to ${placeName(result.place)}: ${e?.message ?? e}`);
    return null;
  }
}
