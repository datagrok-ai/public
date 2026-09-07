/** Node contracts the Bio and Curves catalog entries must honour on a canvas. */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, createNode} from '../rete/node-factory';
import {findNodeTypesProducingOutput} from '../rete/node-factory';
import {PropertyPanel} from '../panel/property-panel';
import {effectiveFuncInputs} from '../utils/func-input-overrides';
import {INCLUDED_FUNC_NQNAMES} from '../rete/included-funcs';
import {inputValueProperty, buildInputValueEditor, resolveInputValue} from '../utils/input-values';
import {emitScript} from '../compiler/script-emitter';
import {hostsHelmEditor, editorBoxSize,
  INLINE_HELM_WIDTH, INLINE_HELM_HEIGHT} from '../rete/scheme';
import {estimateNodeHeight} from '../rete/graph-layout';
import {isLiteralChoiceList} from '../utils/choice-refs';
import {tid} from '../utils/test-ids';
import {makeEditor, destroyEditor, addNode, until} from './test-utils';

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

  test('an editor reporting a rich value still lands the string (the HELM dialog OK path)', async () => {
    // Helm's registered value editor holds a SeqValueBase — its `value` is an object,
    // and the string form lives in `stringValue`. Storing `value` raw serialized as
    // "[object Object]"; combined with the editor never firing changed on the dialog's
    // OK (fixed in Helm), this is why a sketched HELM read back as an empty output.
    const node = createNode('Inputs/Helm Input')!;
    let edits = 0;
    const ed = buildInputValueEditor(node, () => edits++)!;
    try {
      const helm = 'PEPTIDE1{A.C.D}$$$$';
      const input = ed.input!;
      // Helm's editor materializes asynchronously (the platform binds the JS input
      // proxy once the package loads) — poll until the set lands. Programmatic — silent.
      const ready = await until(() => {
        try {
          input.stringValue = helm;
          return true;
        } catch {
          return false;
        }
      }, 20000);
      expect(ready, true, 'the value editor materialized');
      expect(edits, 0, 'a programmatic set is not an edit');
      input.fireChanged(); // exactly what the Helm editor dialog's OK now fires
      expect(await until(() => edits === 1, 3000), true, 'the change was reported once');
      expect(String(node.properties['defaultValue'] ?? ''), helm,
        'the HELM string (never "[object Object]") lands in the configured value');
      const r = resolveInputValue(node);
      expect(r.ok, true, 'the run can feed it to the prepared call');
      expect(r.value, helm);
    } finally {
      ed.root.remove();
    }
  });

  test('the Helm node hosts its editor in a resizable in-card box that it actually fills', async () => {
    const e = makeEditor();
    try {
      const node = await addNode(e.flow, 'Inputs/Helm Input', 100, 100);
      expect(hostsHelmEditor(node), true, 'Macromolecule routes to the helm box');
      const boxSel = `.ff-node[data-node-id="${node.id}"] [data-testid="ff-helm-box"]`;
      expect(await until(() => e.container.querySelector(boxSel) != null, 10000), true,
        'the box renders inside the card — no portal involved');
      expect(e.container.querySelector('.ff-node-preview-portal'), null, 'no portal for value editors');
      const box = e.container.querySelector<HTMLElement>(boxSel)!;
      expect(box.style.width, `${INLINE_HELM_WIDTH}px`, 'the HELM default box');
      expect(box.style.height, `${INLINE_HELM_HEIGHT}px`);
      // A user resize (native CSS handle writes inline style) persists into the
      // node properties and drives the layout estimates.
      box.style.width = '480px';
      box.style.height = '400px';
      expect(await until(() => editorBoxSize(node).width === 480, 5000), true,
        'the resize landed in the node properties');
      expect(estimateNodeHeight(node) > 400, true, 'the card estimate follows the box');
      // THE reported bug: the HELM editor pinned itself to 250×250 and never
      // tracked its container. Once it materializes (async package load), its
      // host must fill the box — at the resized 480px, not the built-in 250.
      const helmMounted = await until(() =>
        e.container.querySelector(`${boxSel} .ui-input-helm`) != null, 20000);
      if (!helmMounted) {
        console.warn('Flow: helm input: the HELM editor did not materialize — fill-the-box unverified');
        return;
      }
      const filled = await until(() => {
        const host = e.container.querySelector<HTMLElement>(`${boxSel} .ui-input-helm .ui-input-editor`);
        return host != null && host.clientWidth >= box.clientWidth - 10;
      }, 10000);
      const host = e.container.querySelector<HTMLElement>(`${boxSel} .ui-input-helm .ui-input-editor`);
      expect(filled, true,
        `the editor host tracks the box (${host?.clientWidth}px vs box ${box.clientWidth}px — was stuck at 250)`);
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 60000});

  test('a value loaded before the editor materializes still renders (async-bind retry)', async () => {
    // The reported bug: a flow loaded with a HELM in `defaultValue` ran fine
    // (the value channel reads the node), but the helm editor itself rendered
    // empty — the initial sync landed on the still-unbound async editor proxy,
    // threw, and nothing retried once the editor materialized.
    const node = createNode('Inputs/Helm Input')!;
    const helm = 'PEPTIDE1{A.C.D}$$$$';
    node.properties['defaultValue'] = helm; // exactly what deserializeFlow leaves behind
    const ed = buildInputValueEditor(node, () => {}, {host: 'node'})!;
    document.body.appendChild(ed.root);
    try {
      const landed = await until(() => {
        try {
          return ed.input!.stringValue === helm;
        } catch {
          return false;
        }
      }, 20000);
      expect(landed, true,
        'the stored value reaches the editor once the async proxy binds — a loaded flow used to show an empty helm');
    } finally {
      ed.root.remove();
    }
  }, {timeout: 30000});

  test('a HELM value runs end to end — braces never ride the header', async () => {
    // The platform's ScriptParser (`_validateParamLine`) takes the span from the
    // line's FIRST `{` to its LAST `}` as the options block, with no string-awareness:
    // `//input: string sequence = "PEPTIDE1{A.C}$$$$" {semType: Macromolecule}` parses
    // its options as `A.C}$$$$" {semType: Macromolecule` and errors out. A
    // brace-carrying default must stay out of the header; every run path goes
    // through DG.Script.create, so this breaks Run/autorun/save alike.
    const e = makeEditor();
    try {
      const helm = 'PEPTIDE1{A.C.D}$$$$';
      const input = await addNode(e.flow, 'Inputs/Helm Input');
      input.properties['defaultValue'] = helm;
      const out = await addNode(e.flow, 'Outputs/Value Output');
      await e.flow.addConnectionByKeys(input.id, 'sequence', out.id, 'value');
      const script = emitScript(e.flow, {name: 'HelmE2E', description: 'test', tags: ['funcflow']});
      const line = script.split('\n').find((l) => l.startsWith('//input:'))!;
      expect(line.includes('{semType: Macromolecule}'), true, 'the qualifier block is intact');
      expect(line.includes('PEPTIDE1'), false, 'the HELM default is omitted, not mangled in place');
      const func = DG.Script.create(script);
      const seq = func.inputs.find((p) => p.name === 'sequence');
      expect(seq != null, true, 'the platform parsed the input parameter');
      expect(seq!.semType, 'Macromolecule', 'the semType survived the parse');
      const fc = func.prepare({sequence: helm});
      await fc.call(undefined, undefined, {processed: true});
      expect(fc.outputs['result'], helm, 'the configured value flowed through the whole run');
    } finally {
      destroyEditor(e);
    }
  }, {timeout: 30000});
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
