import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {readDataframe} from './utils';

/**
 * Regression tests for the canonical FuncCallEditor conversions in Bio.
 *
 * The platform hosts a FuncCallEditor in a dialog and re-evaluates `editor.isValid` on every
 * `onInputChanged` emission. If `isValid` (directly or indirectly) emits `onInputChanged`, the two
 * feed each other and the page freezes. These tests reproduce that platform behavior with a
 * re-entrancy guard and also assert each editor builds (constructor is a common crash site).
 */

/** Simulates the platform's dialog host: on every onInputChanged it re-reads isValid. Returns the
 * maximum re-entrancy depth observed. Bounded so a looping editor terminates the test instead of
 * hanging the whole run. */
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

/** Prepares the target funccall and invokes its editor function through the platform (via
 * grok.functions.call), exactly as the runtime does: target's `editor:` tag → editor func → DG.Widget. */
async function buildEditor(targetName: string, editorName: string, inputs: object): Promise<DG.FuncCallEditor> {
  const call = DG.Func.find({package: 'Bio', name: targetName})[0].prepare(inputs);
  const widget = await grok.functions.call(`Bio:${editorName}`, {call});
  return widget as DG.FuncCallEditor;
}

category('Bio editors: no infinite loop', () => {
  let view: DG.TableView;

  before(async () => {
    // getSeqHelper() runs Bio:getSeqHelper -> initBio through the platform, initializing the
    // package singleton (_package.seqHelper) that GetRegionEditor relies on
    await getSeqHelper();
    const df = await readDataframe('tests/100_3_clustests.csv');
    await grok.data.detectSemanticTypes(df);
    view = grok.shell.addTableView(df);
  });

  after(async () => {
    view?.close();
  });

  test('Sequence Space editor does not loop', async () => {
    const editor = await buildEditor('sequenceSpaceTopMenu', 'SequenceSpaceEditor', {table: view.dataFrame});
    const maxDepth = measureReentrancy(editor, () => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
    });
    expect(maxDepth <= 2, true,
      `Sequence Space editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  test('Seq Activity Cliffs editor does not loop', async () => {
    const editor = await buildEditor('activityCliffs', 'SeqActivityCliffsEditor', {table: view.dataFrame});
    const maxDepth = measureReentrancy(editor, () => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
    });
    expect(maxDepth <= 2, true,
      `Seq Activity Cliffs editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  test('Get Region editor does not loop', async () => {
    const editor = await buildEditor('getRegionTopMenu', 'GetRegionEditor', {table: view.dataFrame});
    const maxDepth = measureReentrancy(editor, () => {
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      const _ = editor.isValid;
    });
    expect(maxDepth <= 2, true,
      `Get Region editor isValid re-entered ${maxDepth} times (expected <=2) — infinite loop`);
  });

  test('Get Region editor history restore', async () => {
    const editor: any = await buildEditor('getRegionTopMenu', 'GetRegionEditor', {table: view.dataFrame});
    // a valid start/end (positions that exist in the sequence column) and a custom name are restored
    editor.loadHistoryString(JSON.stringify({start: '2', end: '5', name: 'my region'}));
    expect(editor.getHistoryString().length > 0, true, 'Get Region editor should expose history');
    const restored = JSON.parse(editor.getHistoryString());
    expect(restored.name, 'my region', 'name not restored');
    expect(restored.start, '2', 'start not restored');
    expect(restored.end, '5', 'end not restored');
    // an out-of-range position must be ignored (kept at the current value)
    const before = JSON.parse(editor.getHistoryString());
    editor.loadHistoryString(JSON.stringify({start: 'no-such-pos', end: 'no-such-pos', name: before.name}));
    const after = JSON.parse(editor.getHistoryString());
    expect(after.start, before.start, 'out-of-range start must be ignored');
    expect(after.end, before.end, 'out-of-range end must be ignored');
  });
});
