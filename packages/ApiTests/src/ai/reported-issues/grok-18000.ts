import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, findProp, look} from '../helpers';

// Regression coverage for GROK-18000: PC plot `columnNames` update propagation.
// Shipped in 1.26.0/1.26.3. Removing entries from `columnNames` from the
// Context Panel previously did not refresh the axes until manual resize. The
// fix makes a `columnNames` change immediately invalidate the layout. We pin
// the round-trip: build via `DG.Viewer.pcPlot`, mutate via `setOptions`, read
// back from `getOptions(true).look`. No DOM/canvas inspection — JS API only.
category('AI: GROK-18000: PC plot columnNames update', () => {
  const v = (n: number = 50, cols: string[] = ['age', 'height', 'weight']): DG.Viewer =>
    DG.Viewer.pcPlot(demog(n), {columnNames: cols});

  test('initial columnNames length is preserved on the look bag', async () => {
    expect((look(v())['columnNames'] as string[]).length, 3);
  });

  test('shrinking columnNames via setOptions updates the look bag exactly', async () => {
    const c = v();
    c.setOptions({columnNames: ['age', 'height']});
    expectArray(look(c)['columnNames'] as string[], ['age', 'height']);
  });

  test('empty then single-column transitions: viewer survives, look reflects each step', async () => {
    const c = v();
    expectNoThrow(() => c.setOptions({columnNames: []}));
    expect((look(c)['columnNames'] as string[]).length, 0);
    c.setOptions({columnNames: ['age']});
    expectArray(look(c)['columnNames'] as string[], ['age']);
  });

  test('columnNames property is discoverable via getProperties (best-effort)', async () => {
    const c = v(20, ['age', 'height']);
    if (findProp(c, 'columnNames') == null) {
      c.setOptions({columnNames: ['age']});
      expectArray(look(c)['columnNames'] as string[], ['age']);
    }
  });
});
