import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';
import {runScaffoldHopping} from '../analysis/scaffold-hopping/scaffold-hopping';

/** BCR-ABL kinase inhibitors — the same 5-row test set the implementation's
 *  composite-weight comment in `scaffold-hopping.ts` cites as its empirical
 *  anchor. Imatinib is the reference; the other four are increasingly distant
 *  analogues:
 *  - Nilotinib: close analog (same diaminopyrimidine + benzamide + methylpiperazine
 *    shape, swapped methylpiperazine for imidazole + trifluoromethyl) → high Tc,
 *    Maeda atom-ratio above 0.4 → expected NOT flagged at Hard preset.
 *  - Dasatinib: thiazole + aminopyrimidine, very different scaffold from Imatinib
 *    → low Tc, Maeda atom-ratio below 0.4 → expected flagged at Hard preset.
 *  - Bosutinib: quinoline-based, also very different → similar to Dasatinib.
 *  - Ponatinib: imidazo[1,2-b]pyridazine, large topological hop → low Tc,
 *    Maeda atom-ratio below 0.4 → expected flagged at Hard preset. */
const BCR_ABL_SMILES = [
  'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1',
  'Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc2C(F)(F)F)cn1',
  'Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1',
  'COc1cc(Nc2ncnc3cc(OCCCN4CCN(C)CC4)c(OC)cc23)c(Cl)cc1Cl',
  'Cc1ccc(C#Cc2cnc3cccc(C(=O)Nc4cc(C(F)(F)F)ccc4N4CCN(C)CC4)n23)cc1',
];
const BCR_ABL_NAMES = ['Imatinib', 'Nilotinib', 'Dasatinib', 'Bosutinib', 'Ponatinib'];

const COL_TANIMOTO = 'Scaffold Hop Tanimoto';
const COL_CATS = 'Scaffold Hop CATS Sim';
const COL_MCS_RATIO = 'Scaffold Hop MCS Ratio';
const COL_SCORE = 'Scaffold Hop Score';
const COL_FLAG = 'Scaffold Hop';

function createBcrAblTable(): {table: DG.DataFrame; molecules: DG.Column} {
  const smilesCol = DG.Column.fromStrings('smiles', BCR_ABL_SMILES);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  const namesCol = DG.Column.fromStrings('name', BCR_ABL_NAMES);
  const table = DG.DataFrame.fromColumns([namesCol, smilesCol]);
  table.name = 'bcr-abl-smoke';
  return {table, molecules: smilesCol};
}

category('scaffold-hopping smoke', () => {
  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  // 1. Smoke — basic flow on the BCR-ABL set adds the five expected columns
  //    and produces non-NaN Tanimoto for every row. This is the lowest-bar
  //    "did the function run end-to-end" check.
  test('columns_added', async () => {
    const {table, molecules} = createBcrAblTable();
    grok.shell.addTableView(table);
    await runScaffoldHopping(
      table, molecules, 0,
      0.05, 0.7, // wide Tc window so all 4 candidates survive
      0.5, 0.0, // permissive Maeda + permissive CATS so the run produces flags
      '[]', // no marked atoms
      false, false, // pure Maeda — flag depends only on atom-ratio
    );
    expect(table.col(COL_TANIMOTO) !== null, true, `${COL_TANIMOTO} column missing`);
    expect(table.col(COL_CATS) !== null, true, `${COL_CATS} column missing`);
    expect(table.col(COL_MCS_RATIO) !== null, true, `${COL_MCS_RATIO} column missing`);
    expect(table.col(COL_SCORE) !== null, true, `${COL_SCORE} column missing`);
    expect(table.col(COL_FLAG) !== null, true, `${COL_FLAG} column missing`);

    const tcCol = table.col(COL_TANIMOTO)!;
    for (let i = 0; i < table.rowCount; i++) {
      const v = tcCol.get(i);
      expect(typeof v === 'number' && !Number.isNaN(v), true,
        `Tanimoto should be numeric at row ${i}, got ${v}`);
    }
  });

  // 2. Reference-row pinning — the reference's own metrics are forced to 1.0
  //    even with parameter combos that would otherwise leave them as NaN.
  //    Guards the early-return path in `scaffold-hopping.ts:127-135`.
  test('reference_pinning', async () => {
    const {table, molecules} = createBcrAblTable();
    grok.shell.addTableView(table);
    await runScaffoldHopping(
      table, molecules, 0,
      0.05, 0.7,
      0.5, 0.8,
      '[]', true, true,
    );
    const tc = table.col(COL_TANIMOTO)!.get(0);
    const cats = table.col(COL_CATS)!.get(0);
    const mcsRatio = table.col(COL_MCS_RATIO)!.get(0);
    const score = table.col(COL_SCORE)!.get(0);
    const isHop = table.col(COL_FLAG)!.get(0);

    expect(tc, 1.0, 'Reference Tanimoto should be 1.0');
    expect(cats, 1.0, 'Reference CATS should be 1.0');
    expect(mcsRatio, 1.0, 'Reference MCS ratio should be 1.0');
    expect(score, 1.0, 'Reference Score should be 1.0');
    expect(isHop, false, 'Reference itself must not be flagged as a hop');
  });

  // 3. No-survivors — impossibly narrow Tc window pre-filters every candidate,
  //    the run completes without throwing, the reference is still pinned,
  //    and no rows are flagged. Guards the early-return at line 137-145.
  test('no_survivors', async () => {
    const {table, molecules} = createBcrAblTable();
    grok.shell.addTableView(table);
    await runScaffoldHopping(
      table, molecules, 0,
      0.99, 1.0, // impossibly tight — only the reference itself can satisfy
      0.4, 0.8,
      '[]', true, true,
    );
    expect(table.col(COL_TANIMOTO) !== null, true, 'Columns must still be added');
    expect(table.col(COL_FLAG)!.get(0), false, 'Reference must not be self-flagged');

    let flaggedCount = 0;
    for (let i = 0; i < table.rowCount; i++)
      if (table.col(COL_FLAG)!.get(i)) flaggedCount++;
    expect(flaggedCount, 0, 'No rows should be flagged when no survivors');

    const refScore = table.col(COL_SCORE)!.get(0);
    expect(refScore, 1.0, 'Reference Score must remain pinned to 1.0');
  });

  // 4. Pure-Maeda mode — `useTcInFlag=false, useCatsInFlag=false`. The flag
  //    must depend only on `MCS_ratio ≤ mcsRatioMax`; raising mcsRatioMax to
  //    1.0 should flag every survivor (the criterion becomes a tautology).
  test('pure_maeda_mode_tautology', async () => {
    const {table, molecules} = createBcrAblTable();
    grok.shell.addTableView(table);
    await runScaffoldHopping(
      table, molecules, 0,
      0.0, 1.0, // accept everything
      1.0, 0.0, // mcsRatio ≤ 1.0 always true; CATS irrelevant
      '[]', false, false, // pure Maeda — only atom-ratio drives flag
    );
    let nonRefFlagged = 0;
    let nonRefMcsComputed = 0;
    for (let i = 1; i < table.rowCount; i++) {
      if (table.col(COL_FLAG)!.get(i)) nonRefFlagged++;
      const m = table.col(COL_MCS_RATIO)!.get(i);
      if (typeof m === 'number' && !Number.isNaN(m)) nonRefMcsComputed++;
    }
    // Every row whose MCS was computed (= every row, since 4 << TOP_N=200)
    // should be flagged when mcsRatioMax=1.0.
    expect(nonRefFlagged, nonRefMcsComputed,
      `Pure-Maeda + mcsRatioMax=1.0: all MCS-computed rows should be flagged ` +
      `(${nonRefFlagged} flagged of ${nonRefMcsComputed} with MCS)`);
  });

  // 5. Hard preset on BCR-ABL — paper-faithful Maeda thresholds. Imatinib
  //    is the reference; we expect at least ONE of the three large-step
  //    analogues (Dasatinib / Bosutinib / Ponatinib) to be flagged because
  //    they share little topology with Imatinib. Nilotinib is the close
  //    analogue and is expected NOT to be flagged (atom-ratio > 0.4).
  test('hard_preset_finds_at_least_one_hop', async () => {
    const {table, molecules} = createBcrAblTable();
    grok.shell.addTableView(table);
    await runScaffoldHopping(
      table, molecules, 0,
      0.05, 0.3, // Hard preset Tc range
      0.4, 0.0, // paper-faithful Maeda; CATS permissive (in case Python script unavailable)
      '[]', true, false, // Tc window in flag; CATS not in flag (graceful re missing Python)
    );
    let hopCount = 0;
    let imatinibFlagged = false;
    for (let i = 0; i < table.rowCount; i++) {
      if (table.col(COL_FLAG)!.get(i)) {
        hopCount++;
        if (i === 0) imatinibFlagged = true;
      }
    }
    expect(imatinibFlagged, false, 'Imatinib (reference) must not be self-flagged');
    expect(hopCount > 0, true,
      `Hard preset must flag at least one of Nilotinib/Dasatinib/Bosutinib/Ponatinib ` +
      `as a scaffold hop of Imatinib; got ${hopCount}`);
  });
});
