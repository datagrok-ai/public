import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {
  detectNonWaterHetatmInstances,
  findPlDiagramColForSource,
  interactionsColForDiagram,
  PL_DIAGRAM_SEM_TYPE,
  PROLIF_INTERACTIONS_TAG,
  PROLIF_SOURCE_TAG,
  renderInteractionBreakdown,
} from '../utils/prolif';


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
  // matching `PL Interactions` column via the `.%prolif-interactions-col`
  // tag set by `runPlBatch`. Tag-based (not name-based) so a future
  // uniquification on column-name collision doesn't break the pairing.
  // Regression here would make the breakdown read from a stale column or
  // null silently.
  // ------------------------------------------------------------------------
  test('interactionsColForDiagram: tag-linked pair resolves correctly', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('PL Diagram', 1),
      DG.Column.string('PL Interactions', 1),
    ]);
    const diagramCol = df.col('PL Diagram')!;
    diagramCol.tags[PROLIF_INTERACTIONS_TAG] = 'PL Interactions';
    expect(interactionsColForDiagram(diagramCol)?.name, 'PL Interactions');
  });

  test('interactionsColForDiagram: returns null when tag missing', async () => {
    // Column has no tag — caller treats this as "not produced by runPlBatch".
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('PL Diagram', 1),
      DG.Column.string('PL Interactions', 1),
    ]);
    expect(interactionsColForDiagram(df.col('PL Diagram')!), null);
  });

  test('interactionsColForDiagram: returns null when tag points at missing column', async () => {
    // Tag exists but the referenced column was removed (manual user delete).
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('PL Diagram', 1),
    ]);
    const diagramCol = df.col('PL Diagram')!;
    diagramCol.tags[PROLIF_INTERACTIONS_TAG] = 'PL Interactions';
    expect(interactionsColForDiagram(diagramCol), null);
  });

  // ------------------------------------------------------------------------
  // findPlDiagramColForSource — drives the inline-diagram cache reuse in the
  // Molecule3D / docking-pose context panels. Without this, clicking a PDB
  // cell after running a batch would silently re-run the Python script (15-30s)
  // instead of pulling the cached HTML from `column.temp`. The match is by
  // (semType === PL_DIAGRAM_SEM_TYPE) AND (.%prolif-source tag === source name)
  // so unrelated rawPng columns from other features don't get claimed.
  // ------------------------------------------------------------------------
  test('findPlDiagramColForSource: matches diagram col by source tag + rawPng semType', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('pdb', 1),
      DG.Column.string('PL Diagram', 1),
    ]);
    const diagramCol = df.col('PL Diagram')!;
    diagramCol.semType = PL_DIAGRAM_SEM_TYPE;
    diagramCol.tags[PROLIF_SOURCE_TAG] = 'pdb';
    expect(findPlDiagramColForSource(df, 'pdb')?.name, 'PL Diagram');
  });

  test('findPlDiagramColForSource: returns null when no batch has run', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.string('pdb', 1)]);
    expect(findPlDiagramColForSource(df, 'pdb'), null);
  });

  test('findPlDiagramColForSource: ignores rawPng cols from unrelated features', async () => {
    // A different rawPng column (e.g. from PowerGrid demo / another package)
    // must not be matched — only ones with the prolif source tag for OUR source.
    const df = DG.DataFrame.fromColumns([
      DG.Column.string('pdb', 1),
      DG.Column.string('other-png', 1),
    ]);
    df.col('other-png')!.semType = PL_DIAGRAM_SEM_TYPE;
    expect(findPlDiagramColForSource(df, 'pdb'), null);
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

});
