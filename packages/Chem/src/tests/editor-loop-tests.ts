import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {createTableView} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';
import {PackageFunctions} from '../package';

/**
 * Regression tests for the canonical FuncCallEditor conversions.
 *
 * The platform hosts a FuncCallEditor in a dialog and re-evaluates `editor.isValid` on every
 * `onInputChanged` emission. If `isValid` (directly or indirectly) emits `onInputChanged`, the
 * two feed each other and the page freezes in an infinite loop. These tests reproduce that
 * platform behavior with a re-entrancy guard: reading `isValid` must not spin.
 */

/** Simulates the platform's dialog host: on every onInputChanged it re-reads isValid.
 * Returns the maximum re-entrancy depth observed while poking the editor. Bounded so a buggy
 * (looping) editor terminates the test instead of hanging the whole run. */
function measureReentrancy(editor: DG.FuncCallEditor, poke: () => void, guard = 200): number {
  let depth = 0;
  let maxDepth = 0;
  const sub = editor.onInputChanged.subscribe(() => {
    depth++;
    maxDepth = Math.max(maxDepth, depth);
    if (depth < guard) {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid; // platform re-validates on each change
    }
    depth--;
  });
  try {
    poke();
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    const _ = editor.isValid; // one more explicit validation, as the platform does before OK
  } finally {
    sub.unsubscribe();
  }
  return maxDepth;
}

/** Prepares the target funccall and invokes its editor function to build the widget, exactly as
 * the platform does (target's `editor:` tag → editor function → DG.Widget). Pass `inputs` to seed
 * the call (e.g. the current table) — editors read `call.inputs['table'] ?? grok.shell.t`, and in
 * the test harness `grok.shell.t` may be null, which would skip the column-input code paths. */
async function buildEditor(targetName: string,
  makeEditor: (call: DG.FuncCall) => DG.Widget | Promise<DG.Widget>,
  inputs: object = {}): Promise<DG.FuncCallEditor> {
  const call = DG.Func.find({package: 'Chem', name: targetName})[0].prepare(inputs);
  return await makeEditor(call) as DG.FuncCallEditor;
}

category('Editors: no infinite loop', () => {
  let view: DG.TableView;

  before(async () => {
    chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
    await chemCommonRdKit.initRdKitModuleLocal();
    view = await createTableView('demo_files/matched_molecular_pairs.csv');
  });

  after(async () => {
    view?.close();
  });

  // The exact repro from the bug report: open the MMP dialog on matched_molecular_pairs.csv and
  // pick the activity columns. Reading isValid must not re-emit onInputChanged.
  test('MMP editor does not loop on activity selection', async () => {
    // eslint-disable-next-line new-cap
    const editor = await buildEditor('mmpAnalysis', (call) => PackageFunctions.MMPEditor(call),
      {table: view.dataFrame});
    const df = view.dataFrame;
    const numericCols = Array.from(df.columns.numerical);
    const maxDepth = measureReentrancy(editor, () => {
      // reading isValid alone must be side-effect free
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
      // emulate the user changing the activities selection, which triggers onInputChanged
      const activitiesInput = (editor as any).activitiesInput as DG.InputBase | undefined;
      if (activitiesInput && numericCols.length)
        activitiesInput.value = numericCols.slice(0, Math.min(2, numericCols.length));
    });
    expect(maxDepth <= 2, true,
      `MMP editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  test('Chem Space editor does not loop', async () => {
    // eslint-disable-next-line new-cap
    const editor = await buildEditor('chemSpaceTopMenu', (call) => PackageFunctions.ChemSpaceEditor(call));
    const maxDepth = measureReentrancy(editor, () => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
    });
    expect(maxDepth <= 2, true,
      `Chem Space editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  test('Activity Cliffs editor does not loop', async () => {
    // eslint-disable-next-line new-cap
    const editor = await buildEditor('activityCliffs', (call) => PackageFunctions.ActivityCliffsEditor(call));
    const maxDepth = measureReentrancy(editor, () => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
    });
    expect(maxDepth <= 2, true,
      `Activity Cliffs editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  test('Map Identifiers editor does not loop', async () => {
    const editor = await buildEditor('getMapIdentifiers',
      (call) => PackageFunctions.mapIdentifiersEditor(call), {table: view.dataFrame});
    const maxDepth = measureReentrancy(editor, () => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
    });
    expect(maxDepth <= 2, true,
      `Map Identifiers editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  test('Descriptors editor does not loop', async () => {
    const editor = await buildEditor('descriptorsDocker',
      (call) => PackageFunctions.descriptorsEditor(call), {table: view.dataFrame});
    const maxDepth = measureReentrancy(editor, () => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
    });
    expect(maxDepth <= 2, true,
      `Descriptors editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  // regression: the constructor used to insertBefore against a not-yet-attached anchor and threw
  // NotFoundError whenever a current table was present
  test('Map Identifiers editor builds with a current table', async () => {
    const editor: any = await buildEditor('getMapIdentifiers',
      (call) => PackageFunctions.mapIdentifiersEditor(call), {table: view.dataFrame});
    expect(editor.columnInput != null, true, 'column input not created');
    expect(editor.root.contains(editor.columnInput.root), true, 'column input not attached to the editor');
    expect(editor.root.contains(editor.fromSourceInput.root), true, 'from-source input not attached');
  });

  test('MMP editor history restore', async () => {
    console.log('[mmp-history] building editor');
    // eslint-disable-next-line new-cap
    const editor: any = await buildEditor('mmpAnalysis', (call) => PackageFunctions.MMPEditor(call),
      {table: view.dataFrame});
    console.log('[mmp-history] editor built');
    const df = view.dataFrame;
    const numCols = Array.from(df.columns.numerical);
    console.log(`[mmp-history] numerical columns: ${numCols.map((c) => c.name).join(', ')}`);
    expect(numCols.length >= 1, true, 'test table has no numerical columns');
    console.log('[mmp-history] loading history string');
    editor.loadHistoryString(JSON.stringify({
      activities: [
        {name: numCols[0].name, diffType: 'ratio', scaling: 'lg'},
        {name: 'no such column', diffType: 'delta', scaling: 'none'},
      ],
      fragmentCutoff: 0.7,
      runOnFilteredData: false,
    }));
    console.log('[mmp-history] history loaded, asserting');
    expect(editor.cutoffInput.value, 0.7, 'cutoff not restored');
    expect(editor.runOnFilteredDataInput.value, false, 'runOnFilteredData not restored');
    console.log('[mmp-history] scalar inputs OK, reading activities');
    const selected: DG.Column[] = editor.activitiesInput.value;
    console.log(`[mmp-history] restored activities: ${selected.map((c) => c.name).join(', ')}`);
    expect(selected.length, 1, 'non-existing column must be dropped on restore');
    expect(selected[0].name, numCols[0].name, 'existing numerical column not restored');
    expect(editor.getDiffTypes()[0], 'ratio', 'diff type not restored');
    const expectedScaling = numCols[0].stats.min < 0 ? 'none' : 'lg';
    expect(editor.getScalings()[0], expectedScaling, 'scaling not restored');
    // the per-activity delta/scaling selector rows must be rebuilt on restore
    const paramRows = (editor.activitiesParamsDiv as HTMLElement).children.length;
    console.log(`[mmp-history] delta/scaling rows rendered: ${paramRows}`);
    expect(paramRows, 1, 'delta/scaling selectors not shown after history restore');
    const roundTrip = JSON.parse(editor.getHistoryString());
    expect(roundTrip.activities.length, 1, 'history round-trip lost the restored activity');
    expect(roundTrip.fragmentCutoff, 0.7, 'history round-trip lost the cutoff');
  });

  test('Map Identifiers editor history restore', async () => {
    const editor: any = await buildEditor('getMapIdentifiers',
      (call) => PackageFunctions.mapIdentifiersEditor(call), {table: view.dataFrame});
    editor.loadHistoryString(JSON.stringify({fromSource: 'chembl', toSource: 'pubchem'}));
    expect(editor.fromSourceInput.value, 'chembl', 'fromSource not restored');
    expect(editor.toSourceInput.value, 'pubchem', 'toSource not restored');
    // a source that is not among the available options must not be applied
    editor.loadHistoryString(JSON.stringify({fromSource: 'no such source', toSource: 'drugbank'}));
    expect(editor.fromSourceInput.value, 'chembl', 'unknown fromSource must be ignored');
    expect(editor.toSourceInput.value, 'drugbank', 'valid toSource not restored');
    const roundTrip = JSON.parse(editor.getHistoryString());
    expect(roundTrip.fromSource, 'chembl', 'history round-trip lost fromSource');
    expect(roundTrip.toSource, 'drugbank', 'history round-trip lost toSource');
  });
// clear: false — the runner closes views between tests by default, and these tests share the
// table view created in before(); editors read grok.shell.tv when the call has no table input
}, {clear: false});
