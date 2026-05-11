// Docking-specific glue around the shared ProLIF panel + batch handler
// (which lives in `@datagrok-libraries/bio/src/prolif/prolif-panel`).
//
// The two Docking-specific things — receptor pre-fetch and the
// `isApplicableAutodock` gate — feed into bio's `runPlBatch` via the
// `prepare` and `buildRowArgs` options. Everything else (UI breakdown,
// capture container, cache, column layout) is shared.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  makeProlifWidget as makeBioProlifWidget,
  runPlBatch,
  ProlifBatchCtx,
} from '@datagrok-libraries/bio/src/prolif/prolif-panel';

import {getReceptorData, RECEPTOR_NAME_RE} from './utils';

export {ProlifBatchCtx};

/** Same heuristic as `Docking:isApplicableAutodock` — duplicated here so
 *  the batch loop can call it without going through the Datagrok function
 *  registry (we'd lose 1-2 RPCs per row). */
function isApplicableAutodock(molecule: string): boolean {
  return molecule.includes('binding energy');
}

/** Docking's variant of the PL batch: pre-fetches the receptor once from the
 *  first valid pose's REMARK header (poses don't carry a full receptor —
 *  ligand and receptor are in separate files paired by REMARK metadata),
 *  then routes the receptor as `protein` and each row's pose as `ligand`.
 *  Skips rows whose pose value doesn't match `isApplicableAutodock`. */
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
      if (!pose || !isApplicableAutodock(pose)) continue;
      if (RECEPTOR_NAME_RE.test(pose)) {
        firstPoseWithReceptor = pose;
        break;
      }
    }
    if (firstPoseWithReceptor == null) {
      grok.shell.warning(
        'Could not resolve receptor for PL batch — no usable poses found.');
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
      if (!pose || (ligandCol != null && !isApplicableAutodock(pose))) return null;
      return {
        protein: receptor ?? pose,
        ligand: ligandCol != null ? pose : '',
        ligand_resname: '',
      };
    },
  });
}

/** Docking variant of `makeProlifWidget` — same UI as bio's default but
 *  wires up the Docking-specific `runDockingPlBatch` runner. */
export function makeProlifWidget(params: {
  protein: string;
  ligand?: string;
  ligand_resname?: string;
}, batchCtx?: ProlifBatchCtx): DG.Widget {
  const runner = batchCtx != null ? () => runDockingPlBatch(batchCtx) : undefined;
  return makeBioProlifWidget(params, batchCtx, runner);
}
