/** Node contracts the Bio and Curves catalog entries must honour on a canvas. */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, createNode} from '../rete/node-factory';
import {findNodeTypesProducingOutput} from '../rete/node-factory';
import {PropertyPanel} from '../panel/property-panel';
import {effectiveFuncInputs} from '../utils/func-input-overrides';
import {INCLUDED_FUNC_NQNAMES} from '../rete/included-funcs';
import {inputValueProperty} from '../utils/input-values';
import {isLiteralChoiceList} from '../utils/choice-refs';
import {tid} from '../utils/test-ids';
import {makeEditor, destroyEditor, addNode} from './test-utils';

function typeNameOf(nqName: string): string | null {
  return getRegisteredFuncs().find((f) => {
    try {
      return f.func.nqName === nqName;
    } catch {
      return false;
    }
  })?.nodeTypeName ?? null;
}

function fastaColumn(name: string, seqs: string[]): DG.Column<string> {
  const col = DG.Column.fromStrings(name, seqs);
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.meta.units = 'fasta';
  col.setTag('aligned', 'SEQ');
  col.setTag('alphabet', 'PT');
  return col;
}

async function expectTableAndColumn(
  nqName: string, columnName: string, semType: string, extra?: string,
): Promise<void> {
  const typeName = typeNameOf(nqName);
  if (!typeName) return; // package not deployed on this stand
  const e = makeEditor();
  try {
    const node = await addNode(e.flow, typeName);
    const inputs = effectiveFuncInputs(node.dgFunc!);
    expect(String(inputs[0].propertyType), 'dataframe', `${nqName}: leads with a table`);
    const col = inputs.find((p) => p.name === columnName);
    expect(col !== undefined, true, `${nqName}: declares ${columnName}`);
    expect(String(col!.propertyType), 'column', `${nqName}: ${columnName} is a real column slot`);
    if (semType)
      expect(col!.semType, semType, `${nqName}: ${columnName} is semType-filtered`);
    const tables = node.properties['columnTables'] as Record<string, string> | undefined;
    expect(tables?.[columnName], inputs[0].name, `${nqName}: ${columnName} resolves against the table`);
    expect(node.requiredInputs.includes(inputs[0].name), true, `${nqName}: table required`);
    expect(node.requiredInputs.includes(columnName), true, `${nqName}: ${columnName} required`);
    if (extra)
      expect(extra in node.inputs, true, `${nqName}: ${extra} slot exists`);
  } finally {
    destroyEditor(e);
  }
}

category('Flow: bio nodes', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('the void / same-frame bio entries are replaced by column-returning twins', async () => {
    const replaced: Record<string, string> = {
      'Bio:getRegionTopMenu': 'Bio:extractRegion',
      'Bio:toAtomicLevel': 'Bio:toAtomicLevelColumn',
      'Bio:moleculesToHelmTopMenu': 'Bio:moleculesToHelmColumn',
      'Bio:splitToMonomersTopMenu': 'Bio:splitToMonomersColumns',
      'Bio:immunumAntibodyNumbering': 'Bio:applyAntibodyNumbering',
    };
    for (const [gone, replacement] of Object.entries(replaced)) {
      expect(INCLUDED_FUNC_NQNAMES.has(gone), false, `${gone} superseded`);
      expect(INCLUDED_FUNC_NQNAMES.has(replacement), true, `${replacement} is the replacement`);
    }
  });

  test('the new bio nodes are on the allowlist', async () => {
    for (const added of ['Bio:convertNotation', 'Bio:tagAsMacromolecule', 'Bio:motifSearch',
      'Bio:GenerateSequences', 'Bio:toAtomicLevelSingleSeq'])
      expect(INCLUDED_FUNC_NQNAMES.has(added), true, `${added} listed`);
    expect(INCLUDED_FUNC_NQNAMES.has('Bio:seq2atomic'), false, 'the duplicate name stays out');
  });

  test('every bio twin leads with a table its sequence column resolves against', async () => {
    await expectTableAndColumn('Bio:extractRegion', 'sequence', 'Macromolecule', 'start');
    await expectTableAndColumn('Bio:toAtomicLevelColumn', 'sequence', 'Macromolecule');
    await expectTableAndColumn('Bio:splitToMonomersColumns', 'sequence', 'Macromolecule');
    await expectTableAndColumn('Bio:convertNotation', 'sequence', 'Macromolecule', 'targetNotation');
    await expectTableAndColumn('Bio:motifSearch', 'sequence', 'Macromolecule', 'motif');
    await expectTableAndColumn('Bio:applyAntibodyNumbering', 'sequence', 'Macromolecule', 'scheme');
    await expectTableAndColumn('Bio:moleculesToHelmColumn', 'molecules', 'Molecule');
  });

  test('Tag as Macromolecule does not filter its input by semType', async () => {
    const typeName = typeNameOf('Bio:tagAsMacromolecule');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const inputs = effectiveFuncInputs(node.dgFunc!);
      const col = inputs.find((p) => p.name === 'sequence')!;
      expect(String(col.propertyType), 'column');
      expect(!col.semType, true, 'no semType filter — a derived column has none yet');
      const tables = node.properties['columnTables'] as Record<string, string> | undefined;
      expect(tables?.['sequence'], inputs[0].name, 'still resolves against the table');
    } finally {
      destroyEditor(e);
    }
  });

  test('the bio twins declare literal choices, not free text', async () => {
    for (const [nqName, param, expected] of [
      ['Bio:convertNotation', 'targetNotation', ['fasta', 'separator', 'helm']],
      ['Bio:applyAntibodyNumbering', 'scheme', ['imgt', 'kabat']],
    ] as [string, string, string[]][]) {
      const typeName = typeNameOf(nqName);
      if (!typeName) continue;
      const func = getRegisteredFuncs().find((f) => f.nodeTypeName === typeName)!.func;
      const prop = effectiveFuncInputs(func).find((p) => p.name === param)!;
      const choices = prop.choices ?? [];
      expect(isLiteralChoiceList(choices), true, `${nqName}.${param} has real choices`);
      for (const v of expected)
        expect(choices.includes(v), true, `${nqName}.${param} offers ${v}`);
    }
  });

  test('the bio twin panel offers a column picker for its sequence column', async () => {
    const typeName = typeNameOf('Bio:convertNotation');
    if (!typeName) return;
    const e = makeEditor();
    const panel = new PropertyPanel(e.flow);
    panel.onPickColumns = () => {}; // the icon is gated on a handler being wired
    document.body.appendChild(panel.root);
    try {
      const node = await addNode(e.flow, typeName);
      panel.showNode(node);
      expect(!!panel.root.querySelector('[data-param="sequence"]'), true, 'the column row renders');
      expect(!!panel.root.querySelector(`[data-testid="${tid('prop-pick-columns', 'sequence')}"]`), true,
        'and carries the picker icon');
      expect(!!panel.root.querySelector('[data-param="targetNotation"]'), true, 'the notation row renders');
    } finally {
      panel.root.remove();
      destroyEditor(e);
    }
  });
});

category('Flow: helm input', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('Helm Input is a string input tagged Macromolecule', async () => {
    const node = createNode('Inputs/Helm Input')!;
    expect(!!node, true, 'the node type is registered');
    expect(node.dgOutputType, 'string', 'a HELM sequence IS a string');
    expect(node.properties['semType'], 'Macromolecule');
    expect('sequence' in node.outputs, true, 'one named output socket');
  });

  test('the built property carries the semType through to the value editor', async () => {
    const node = createNode('Inputs/Helm Input')!;
    const prop = inputValueProperty(node);
    expect(prop !== null, true, 'a string input has an inline value editor');
    expect(prop!.semType, 'Macromolecule', 'without this the editor is a text box');
    expect(String(prop!.propertyType), 'string');
  });

  test('a bare string slot still leads with String Input', async () => {
    const matches = findNodeTypesProducingOutput('string');
    const inputs = matches.filter((m) => m.typeName.startsWith('Inputs/'));
    expect(inputs.length > 0, true, 'string inputs are offered');
    expect(inputs[0].typeName, 'Inputs/String Input', 'the specializations rank below the general one');
    expect(matches.some((m) => m.typeName === 'Inputs/Helm Input'), true, 'Helm Input is still offered');
  });
});

category('Flow: bio functions', () => {
  test('Convert Sequence Notation returns the converted column (needs Bio deployed)', async () => {
    const func = DG.Func.find({package: 'Bio', name: 'convertNotation'})[0];
    if (!func) return;
    const sequence = fastaColumn('seq', ['MDYKETL', 'MDYKQTL']);
    const table = DG.DataFrame.fromColumns([sequence]);
    const col: DG.Column = await func.apply({table, sequence, targetNotation: 'separator', separator: '-'});
    expect(col instanceof DG.Column, true, 'returns a column');
    expect(col.length, 2, 'one converted sequence per row');
    expect(table.col(col.name) !== null, true, 'and it is in the table, so the pass-through carries it');
    expect(String(col.get(0)).includes('-'), true, 'converted to the requested separator notation');
  });

  test('Motif Search filters rows by motif (needs Bio deployed)', async () => {
    const func = DG.Func.find({package: 'Bio', name: 'motifSearch'})[0];
    if (!func) return;
    const sequence = fastaColumn('seq', ['MDYKETL', 'AAAAAAA', 'MDYKQTL']);
    const table = DG.DataFrame.fromColumns(
      [sequence, DG.Column.fromInt32Array('id', new Int32Array([1, 2, 3]))]);
    const res: DG.DataFrame = await func.apply({table, sequence, motif: 'MDYK'});
    expect(res instanceof DG.DataFrame, true, 'returns a dataframe');
    expect(res.rowCount, 2, 'the two MDYK sequences match');
    expect(res.columns.names().includes('id'), true, 'the other columns come along');
    expect(table.rowCount, 3, 'the input table is not mutated');
  });

  test('Tag as Macromolecule tags a plain string column (needs Bio deployed)', async () => {
    const func = DG.Func.find({package: 'Bio', name: 'tagAsMacromolecule'})[0];
    if (!func) return;
    // Deliberately untagged — this is what a column built inside a flow looks like.
    const sequence = DG.Column.fromStrings('peptide', [
      'MDYKETLLMPKTDFPMRGGLPNKEPQI', 'MDYKQTLLMPKTDFPMRGGLPNKEPQI', 'MDYKETLLMPKTDFPMRGGLPNKEPKI']);
    const table = DG.DataFrame.fromColumns([sequence]);
    expect(sequence.semType !== DG.SEMTYPE.MACROMOLECULE, true, 'starts untagged');
    const col: DG.Column = await func.apply({table, sequence});
    expect(col.semType, DG.SEMTYPE.MACROMOLECULE, 'detection ran and tagged it');
    expect(col.name, 'peptide', 'the same column, tagged in place — not a copy');
  });

  test('Split to Monomers returns the columns it created (needs Bio deployed)', async () => {
    const func = DG.Func.find({package: 'Bio', name: 'splitToMonomersColumns'})[0];
    if (!func) return;
    const sequence = fastaColumn('seq', ['MDYK', 'MDYQ']);
    sequence.setTag('aligned', 'SEQ.MSA');
    const table = DG.DataFrame.fromColumns([sequence]);
    const res: DG.DataFrame = await func.apply({table, sequence});
    expect(res instanceof DG.DataFrame, true, 'returns a dataframe');
    expect(res !== table, true, 'a frame of the monomer columns, NOT the mutated input handed back');
    expect(res.columns.length > 0, true, 'one column per position');
    expect(res.rowCount, table.rowCount, 'aligned row-for-row with the input');
    for (const name of res.columns.names())
      expect(table.col(name) !== null, true, `${name} is also appended to the input table`);
  });
});

category('Flow: curves nodes', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('the curve statistic twins are listed, the name-addressed originals are not', async () => {
    for (const added of ['Curves:addCurveStatistic', 'Curves:addAggrCurveStatistic', 'Curves:dataToCurves'])
      expect(INCLUDED_FUNC_NQNAMES.has(added), true, `${added} listed`);
    // Still registered (the Fit pane and Data to Curves call them by name), just not offered as nodes.
    for (const gone of ['Curves:addStatisticsColumn', 'Curves:addAggrStatisticsColumn'])
      expect(INCLUDED_FUNC_NQNAMES.has(gone), false, `${gone} superseded by its column-slot twin`);
  });

  test('the curve twins take a real fit column slot', async () => {
    await expectTableAndColumn('Curves:addCurveStatistic', 'curves', 'fit', 'statistic');
    await expectTableAndColumn('Curves:addAggrCurveStatistic', 'curves', 'fit', 'aggregation');
  });

  test('the free-text statistic and aggregation fields became choices', async () => {
    for (const [nqName, param, expected] of [
      ['Curves:addCurveStatistic', 'statistic', ['rSquared', 'auc', 'interceptX']],
      ['Curves:addAggrCurveStatistic', 'statistic', ['rSquared', 'auc', 'interceptX']],
      ['Curves:addAggrCurveStatistic', 'aggregation', ['med', 'avg', 'min', 'max']],
    ] as [string, string, string[]][]) {
      const typeName = typeNameOf(nqName);
      if (!typeName) continue;
      const func = getRegisteredFuncs().find((f) => f.nodeTypeName === typeName)!.func;
      const prop = effectiveFuncInputs(func).find((p) => p.name === param)!;
      const choices = prop.choices ?? [];
      expect(isLiteralChoiceList(choices), true, `${nqName}.${param} has real choices`);
      for (const v of expected)
        expect(choices.includes(v), true, `${nqName}.${param} offers ${v}`);
    }
  });

  test('Data to Curves constrains its column pickers by type', async () => {
    const typeName = typeNameOf('Curves:dataToCurves');
    if (!typeName) return;
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, typeName);
      const inputs = effectiveFuncInputs(node.dgFunc!);
      const tables = node.properties['columnTables'] as Record<string, string> | undefined;
      for (const name of ['concentrationCol', 'readoutCol', 'batchIDCol', 'assayCol',
        'runIDCol', 'compoundIDCol', 'targetEntityCol']) {
        const p = inputs.find((i) => i.name === name);
        expect(p !== undefined, true, `declares ${name}`);
        expect(tables?.[name], inputs[0].name, `${name} resolves against the table`);
        expect(node.requiredInputs.includes(name), true, `${name} required`);
      }
      expect(node.requiredInputs.includes('parentTable'), false,
        'the parent-level inputs are optional — the node runs on well-level data alone');
    } finally {
      destroyEditor(e);
    }
  });

  test('Add Curve Statistic extracts a statistic (needs Curves deployed)', async () => {
    const func = DG.Func.find({package: 'Curves', name: 'addCurveStatistic'})[0];
    if (!func) return;
    // Probe the sibling first: a stand can serve a bundle older than the registered metadata,
    // and the sibling run distinguishes that from this function being broken.
    const sibling = DG.Func.find({package: 'Curves', name: 'addStatisticsColumn'})[0];
    const curve = JSON.stringify({
      series: [{fitFunction: 'linear', points: [{x: 1, y: 2}, {x: 2, y: 4}, {x: 3, y: 6}, {x: 4, y: 8}]}],
    });
    const curves = DG.Column.fromStrings('curve', [curve, curve]);
    curves.semType = 'fit';
    const table = DG.DataFrame.fromColumns([curves]);
    let siblingRuns = false;
    if (sibling) {
      try {
        const probe = DG.DataFrame.fromColumns([curves.clone()]);
        await sibling.apply({table: probe, colName: 'curve', propName: 'rSquared', seriesNumber: 0});
        siblingRuns = true;
      } catch {/* the whole package is unavailable — the assertions below say so */}
    }

    let col: DG.Column;
    try {
      col = await func.apply({table, curves, statistic: 'rSquared', seriesNumber: 0});
    } catch (e) {
      // Skip only when the package executes but the served bundle lacks this export
      // (a stand pinned to an older Curves build) — asserting further would test the stand.
      const stale = siblingRuns && String((e as Error)?.message ?? e).includes('is not a function');
      if (stale) {
        console.warn('Flow: skipping Add Curve Statistic — this stand serves a Curves bundle ' +
          'older than the registered metadata. Republish Curves as the current version to cover it.');
        return;
      }
      throw e;
    }
    expect(col instanceof DG.Column, true, 'returns a column');
    expect(col.length, 2, 'one value per row');
    expect(table.col(col.name) !== null, true, 'and it is in the table');
    expect(col.get(0) > 0.99, true, 'rSquared of an exact linear fit is ~1');
  });
});
