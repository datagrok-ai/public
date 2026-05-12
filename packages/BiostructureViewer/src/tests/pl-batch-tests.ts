import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {
  detectNonWaterHetatmInstances,
  getPlHtmlForRow,
  interactionsColForDiagram,
  renderInteractionBreakdown,
  runPlBatch,
  PL_DIAGRAM_SEM_TYPE,
  type ProlifBatchCtx,
} from '../utils/prolif';


// Minimal valid PDB used by the end-to-end batch smoke test below.
// Real PDB needs both ATOM (protein backbone) and a non-water HETATM
// (small ligand). The script's binding-site extraction picks the
// neighbourhood around the HETATM — coordinates are colocated here so
// the radius cut doesn't drop everything.
const MINIMAL_PDB =
  'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n' +
  'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C  \n' +
  'ATOM      3  C   ALA A   1       2.000   1.420   0.000  1.00  0.00           C  \n' +
  'HETATM    4  C1  LIG A 100       3.000   1.420   0.000  1.00  0.00           C  \n' +
  'END';


category('PLBatch', () => {
  // ------------------------------------------------------------------------
  // hasNonWaterHetatm — precondition guard for the per-row PL batch loop.
  // Locks in the contract that the batch handler in package.ts uses to
  // decide which rows to process.
  // ------------------------------------------------------------------------
  test('hasNonWaterHetatm: ligand HETATM passes', async () => {
    const pdb =
      'ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N  \n' +
      'HETATM    1  C1  LIG A 100      0.000   0.000   0.000  1.00  0.00           C  ';
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: pdb}) as boolean;
    expect(r, true);
  });

  test('hasNonWaterHetatm: water-only HETATM is rejected', async () => {
    const pdb =
      'ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N  \n' +
      'HETATM    1  O   HOH A 200      0.000   0.000   0.000  1.00  0.00           O  ';
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: pdb}) as boolean;
    expect(r, false);
  });

  test('hasNonWaterHetatm: empty input is rejected', async () => {
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: ''}) as boolean;
    expect(r, false);
  });

  test('hasNonWaterHetatm: AutoDock-pose marker is rejected', async () => {
    // The `binding energy` substring check fires before any other guard,
    // so this returns false even if the rest of the input would otherwise
    // qualify. ATOM/HETATM lines are present here to make the input
    // realistic; the assertion exercises the marker short-circuit.
    const pose =
      'ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N  \n' +
      'REMARK  Name = test\n' +
      'REMARK  binding energy = -5.0\n' +
      'HETATM    1  C1  LIG A 100      0.000   0.000   0.000  1.00  0.00           C  ';
    const r = await grok.functions.call(
      `${_package.name}:hasNonWaterHetatm`, {molecule: pose}) as boolean;
    expect(r, false);
  });

  // ------------------------------------------------------------------------
  // detectNonWaterHetatmInstances — drives the ligand-picker / summary.
  // Two contracts: (1) cofactors filtered out (HEM/NAD/... would otherwise
  // outvote a small inhibitor on the most-common-HETATM heuristic);
  // (2) same resname at multiple chain/resid sites yields multiple entries.
  // ------------------------------------------------------------------------
  test('detectNonWaterHetatmInstances: filters cofactors, lists ligand instances', async () => {
    // PDB columns: 13-16 atom, 17-20 resname, 21 chain, 22-26 resid (1-indexed).
    const pdb =
      'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n' +
      'HETATM    1  FE  HEM A 500       0.000   0.000   0.000  1.00  0.00          FE  \n' +
      'HETATM    2  C1  HEM A 500       1.000   0.000   0.000  1.00  0.00           C  \n' +
      'HETATM    3  C1  STI A 401       2.000   0.000   0.000  1.00  0.00           C  \n' +
      'HETATM    4  C1  STI B 402       3.000   0.000   0.000  1.00  0.00           C  ';
    const result = detectNonWaterHetatmInstances(pdb);
    expect(result.length, 2);
    expect(result.includes('STI A 401'), true);
    expect(result.includes('STI B 402'), true);
    expect(result.some((s) => s.startsWith('HEM')), false);
  });

  // ------------------------------------------------------------------------
  // interactionsColForDiagram — wires the post-batch context panel to the
  // matching `PL Interactions` column. New runs always produce canonical names,
  // but the function does suffix-aware lookup so DataFrames carrying `(N)`-
  // suffixed columns from older runs keep working. Regression here would make
  // the breakdown read from a stale column or null silently.
  // ------------------------------------------------------------------------
  test('interactionsColForDiagram: matches bare PL Diagram column', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('PL Diagram', 1),
      DG.Column.string('PL Interactions', 1),
    ]);
    const col = interactionsColForDiagram(df, 'PL Diagram');
    expect(col?.name, 'PL Interactions');
  });

  test('interactionsColForDiagram: matches suffixed PL Diagram column', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('PL Diagram (2)', 1),
      DG.Column.string('PL Interactions (2)', 1),
    ]);
    const col = interactionsColForDiagram(df, 'PL Diagram (2)');
    expect(col?.name, 'PL Interactions (2)');
  });

  test('interactionsColForDiagram: returns null for non-matching name', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('something else', 1),
    ]);
    expect(interactionsColForDiagram(df, 'something else'), null);
  });

  test('interactionsColForDiagram: returns null when suffix has no sibling', async () => {
    // PL Diagram exists but matching PL Interactions does not — the batch
    // would never produce this state, but a manual column delete could.
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('PL Diagram', 1),
    ]);
    expect(interactionsColForDiagram(df, 'PL Diagram'), null);
  });

  // ------------------------------------------------------------------------
  // renderInteractionBreakdown — parses the comma-separated `RESNAME+RESID_CODE`
  // string the Python script emits and groups residues by interaction type
  // for the post-batch context panel. Three contracts: (1) empty input
  // shows the empty-state message; (2) known codes get human-readable
  // labels and residues sort by numeric id within a code; (3) unknown
  // codes surface under their raw name rather than being silently dropped.
  // ------------------------------------------------------------------------
  test('renderInteractionBreakdown: empty string shows empty-state message', async () => {
    const host = renderInteractionBreakdown('');
    expect(host.textContent!.includes('No interactions detected'), true);
    expect(host.textContent!.includes('(0)'), true);
  });

  test('renderInteractionBreakdown: groups by code, sorts residues by numeric id', async () => {
    // Mixed input: 3 hydrophobic residues out of numeric order + 2 H-bond
    // donors. Locks in (a) labels use the human-readable form ("Hydrophobic"
    // not "HY") and (b) numeric residue sort beats alphabetic (GLY3 before
    // GLN21 before TYR123 — plain string sort would give GLN21 first).
    const host = renderInteractionBreakdown(
      'TYR123_HY, GLN21_HY, GLY3_HY, ASP40_HD, GLU50_HD');
    const text = host.textContent!;
    expect(text.includes('Hydrophobic'), true);
    expect(text.includes('H-bond donors'), true);
    const gly3 = text.indexOf('GLY3');
    const gln21 = text.indexOf('GLN21');
    const tyr123 = text.indexOf('TYR123');
    expect(gly3 > 0 && gly3 < gln21 && gln21 < tyr123, true);
  });

  test('renderInteractionBreakdown: unknown code surfaces under its raw name', async () => {
    // Defensive: if Python adds a new interaction type before TS knows
    // about it, the residue must still show up — under the raw code —
    // rather than being silently dropped.
    const host = renderInteractionBreakdown('ASP40_NOVEL');
    const text = host.textContent!;
    expect(text.includes('NOVEL'), true);
    expect(text.includes('ASP40'), true);
  });

  // ------------------------------------------------------------------------
  // runPlBatch end-to-end smoke. Calls the actual Python script through
  // the Datagrok function registry on a single-row DataFrame, then asserts
  // that the always-present output columns exist with the right semType and
  // that the HTML cache was populated. Marked skipReason because the script
  // worker (conda env with ProLIF, MDAnalysis, RDKit, pdbfixer, OpenBabel)
  // isn't available in every test environment — clear the skipReason locally
  // to actually run it.
  //
  // Note: the Python script's internal helpers (`is_pdbqt_lines`,
  // `pdbqt_to_pdb_lines`, `merge_protein_and_ligand_lines`, etc.) are
  // exercised transitively by this test on a representative PDB. Direct
  // unit tests of those helpers would require extracting them into a
  // separate importable module — deferred until that refactor lands.
  // ------------------------------------------------------------------------
  test('runPlBatch: end-to-end on minimal PDB populates diagram + cache', async () => {
    const pdbCol = DG.Column.fromStrings('protein', [MINIMAL_PDB]);
    const df = DG.DataFrame.fromColumns([pdbCol]);
    const ctx: ProlifBatchCtx = {df, pdbCol};

    await runPlBatch({
      ctx,
      buildRowArgs: () => ({protein: MINIMAL_PDB, ligand: '', ligand_resname: ''}),
    });

    const diagramCol = df.col('PL Diagram');
    const interactionsCol = df.col('PL Interactions');
    expect(diagramCol != null, true);
    expect(interactionsCol != null, true);
    expect(diagramCol?.semType, PL_DIAGRAM_SEM_TYPE);
    // Cache populated by `_setHtmlForRow` during the batch — opens the
    // post-batch panel without re-running the script.
    expect((getPlHtmlForRow(df, 0) ?? '').length > 0, true);
  }, {skipReason: 'Requires the Python script-worker env (conda + ProLIF + RDKit + MDAnalysis + pdbfixer)'});
});
