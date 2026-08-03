import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, expectNoThrow, findProp, look} from '../helpers';

// Regression coverage for GROK-18550: PC plot normalizeEachColumn stability with log columns set.
category('AI: GROK-18550: PC plot normalizeEachColumn stability', () => {
  const v = (cols: string[] = ['age', 'height', 'weight']): DG.Viewer =>
    DG.Viewer.pcPlot(demog(), {columnNames: cols, logColumnsColumnNames: ['age']});

  test('toggle normalizeEachColumn true→false→true with log columns set: no throw', async () => {
    const c = v();
    expectNoThrow(() => {
      c.setOptions({normalizeEachColumn: false});
      c.setOptions({normalizeEachColumn: true});
    });
    expectLook(c, {normalizeEachColumn: true});
    const logs = look(c)['logColumnsColumnNames'] as string[] | undefined;
    if (logs !== undefined)
      expectArray(logs, ['age']);
  });

  test('with empty logColumnsColumnNames: toggle normalizeEachColumn survives', async () => {
    const c = v();
    expectNoThrow(() => {
      c.setOptions({logColumnsColumnNames: []});
      c.setOptions({normalizeEachColumn: false});
      c.setOptions({normalizeEachColumn: true});
      c.setOptions({normalizeEachColumn: false});
    });
    expectLook(c, {normalizeEachColumn: false});
    const logs = look(c)['logColumnsColumnNames'] as string[] | undefined;
    if (logs !== undefined)
      expect(logs.length, 0);
  });

  test('getProperties lists normalizeEachColumn and logColumnsColumnNames', async () => {
    const c = DG.Viewer.pcPlot(demog(20), {columnNames: ['age', 'height'], logColumnsColumnNames: ['age']});
    expect(findProp(c, 'normalizeEachColumn') != null, true);
    expect(findProp(c, 'logColumnsColumnNames') != null, true);
  });
});
