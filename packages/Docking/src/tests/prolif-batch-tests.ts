import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {runPlBatch} from '@datagrok-libraries/bio/src/prolif/prolif-panel';


// Docking's makeProlifWidget wires a `prepare` callback into runPlBatch that
// pre-fetches a single receptor from the first usable AutoDock pose. When no
// usable pose is found (or the file fetch fails), `prepare` returns false and
// the batch must:
//   * skip every per-row script call (no Python work)
//   * remove the two output columns (PL Diagram, PL Interactions) it added
//     up front, so the DataFrame is left in its pre-batch state
//
// If this contract regresses, the user sees empty diagram/interactions
// columns appear and stay — silent partial output.
category('ProlifBatch', () => {
  test('runPlBatch: prepare returning false drops output columns and skips rows', async () => {
    const pdbCol = DG.Column.fromStrings('pose', ['REMARK pose A', 'REMARK pose B']);
    const df = DG.DataFrame.fromColumns([pdbCol]);

    let buildRowArgsCalls = 0;
    await runPlBatch({
      ctx: {df, pdbCol, ligandCol: pdbCol},
      prepare: async () => false,
      buildRowArgs: () => {
        buildRowArgsCalls++;
        return {protein: '', ligand: '', ligand_resname: ''};
      },
    });

    // Per the contract — columns were added up front (so the grid would show
    // them ticking through), then removed when prepare bailed.
    expect(df.col('PL Diagram'), null);
    expect(df.col('PL Interactions'), null);
    // No per-row work happened — buildRowArgs was never invoked.
    expect(buildRowArgsCalls, 0);
  });

  test('runPlBatch: two consecutive prepare-false runs leave only one set of (empty) PL columns', async () => {
    // Re-clicking "Compute for whole dataset" must overwrite the previous
    // PL columns in place rather than append `(2)` / `(3)` suffixed siblings.
    // The bail-via-prepare-false path is the cheapest way to exercise the
    // pre-batch column-drop without standing up a Python script worker.
    const pdbCol = DG.Column.fromStrings('pose', ['REMARK pose A']);
    const df = DG.DataFrame.fromColumns([pdbCol]);
    const noop = async () => false;
    const argsBuilder = () =>
      ({protein: '', ligand: '', ligand_resname: ''});

    // First run — adds PL Diagram + PL Interactions then drops them (prepare
    // returns false). After return, columns are gone.
    await runPlBatch({ctx: {df, pdbCol, ligandCol: pdbCol}, prepare: noop, buildRowArgs: argsBuilder});
    // Manually re-add them to simulate a successful prior run that left
    // canonical columns in the DataFrame, so the next run has something
    // to clean up.
    df.columns.addNewString('PL Diagram');
    df.columns.addNewString('PL Interactions');
    df.columns.addNewInt('PL #H-bond donor');
    const colCountBefore = df.columns.length;

    // Second run — should drop the prior canonical PL columns first, then
    // bail (prepare returns false), leaving the DataFrame in the same shape
    // as before the prior run.
    await runPlBatch({ctx: {df, pdbCol, ligandCol: pdbCol}, prepare: noop, buildRowArgs: argsBuilder});
    expect(df.col('PL Diagram'), null);
    expect(df.col('PL Interactions'), null);
    expect(df.col('PL #H-bond donor'), null);
    expect(df.col('PL Diagram (2)'), null);
    expect(df.col('PL Interactions (2)'), null);
    // Net column count drops by the 3 prior PL columns we manually added.
    expect(df.columns.length, colCountBefore - 3);
  });

  test('runPlBatch: prepare throwing is treated the same as returning false', async () => {
    // Docking's prepare wraps `getReceptorData` in try/catch and routes a
    // network error to a warning + return false. The shared batch handler
    // also catches an unhandled throw from prepare as a safety net — this
    // test locks that net in so a future refactor can't accidentally let
    // exceptions escape and leave the half-added columns in the DataFrame.
    const pdbCol = DG.Column.fromStrings('pose', ['REMARK pose A']);
    const df = DG.DataFrame.fromColumns([pdbCol]);

    let buildRowArgsCalls = 0;
    await runPlBatch({
      ctx: {df, pdbCol, ligandCol: pdbCol},
      prepare: async () => { throw new Error('simulated receptor fetch failure'); },
      buildRowArgs: () => {
        buildRowArgsCalls++;
        return {protein: '', ligand: '', ligand_resname: ''};
      },
    });

    expect(df.col('PL Diagram'), null);
    expect(df.col('PL Interactions'), null);
    expect(buildRowArgsCalls, 0);
  });
});
