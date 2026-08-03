// AutoDock-pose flavour of the ProLIF widget + batch.
//
// The AutoDock-specific bits (everything else is shared via the generic
// batch handler in `./prolif`):
//   * Pre-fetches the receptor from `System:AppData/Docking/targets/`
//     via `getReceptorData` — poses don't carry the receptor themselves.
//   * Gates per-row work by `isAutoDockPose` (the `binding energy`
//     REMARK heuristic).

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  isAutoDockPose,
  makeProlifWidget,
  runPlBatch,
  ProlifBatchCtx,
  ProlifWidgetOptions,
} from './prolif';

import {getReceptorData, RECEPTOR_NAME_RE} from './autodock-receptor';

/** AutoDock-pose batch: pre-fetches the receptor once from the first
 *  valid pose's REMARK header (poses don't carry a full receptor — ligand
 *  and receptor are in separate files paired by REMARK metadata), then
 *  routes the receptor as `protein` and each row's pose as `ligand`.
 *  Skips rows whose pose value doesn't match `isAutoDockPose`. */
async function runDockingPlBatch(ctx: ProlifBatchCtx): Promise<void> {
  const {df, pdbCol, ligandCol} = ctx;
  const rowCount = df.rowCount;

  let receptor: string | null = null;
  const prepare = async (): Promise<boolean> => {
    if (ligandCol == null) return true; // BSV-style fallback — unlikely here
    // Pre-screen poses with a regex match (cheap, no I/O) before the file
    // fetch — getReceptorData makes 3 RPCs per call so we look for the
    // first usable pose first.
    let firstPoseWithReceptor: string | null = null;
    for (let i = 0; i < rowCount; i++) {
      const pose = pdbCol.get(i);
      if (!pose || !isAutoDockPose(pose)) continue;
      if (RECEPTOR_NAME_RE.test(pose)) { firstPoseWithReceptor = pose; break; }
    }
    if (firstPoseWithReceptor == null) {
      grok.shell.warning('Could not resolve receptor for PL batch — no usable poses found.');
      return false;
    }
    try {
      const data = await getReceptorData(firstPoseWithReceptor);
      receptor = typeof data.data === 'string'
        ? data.data
        : new TextDecoder().decode(data.data as Uint8Array);
      return true;
    } catch (err) {
      grok.shell.warning(`Could not fetch receptor for PL batch: ${
        err instanceof Error ? err.message : String(err)}`);
      return false;
    }
  };

  await runPlBatch({
    ctx,
    prepare,
    buildRowArgs: (i) => {
      const pose = pdbCol.get(i);
      if (!pose || (ligandCol != null && !isAutoDockPose(pose))) return null;
      return {
        protein: receptor ?? pose,
        ligand: ligandCol != null ? pose : '',
        ligand_resname: '',
      };
    },
  });
}

/** Docking-flavour wrapper around BSV's `makeProlifWidget`. Wires the
 *  Docking-specific `runDockingPlBatch` runner; extra options pass
 *  through (e.g. `showInlineDiagram` to mount the LigNetwork iframe). */
export function makeDockingProlifWidget(params: {
  protein: string;
  ligand?: string;
  ligand_resname?: string;
}, batchCtx?: ProlifBatchCtx, opts?: ProlifWidgetOptions): DG.Widget {
  const runBatch = batchCtx != null ? () => runDockingPlBatch(batchCtx) : undefined;
  return makeProlifWidget(params, batchCtx, {...opts, runBatch});
}
