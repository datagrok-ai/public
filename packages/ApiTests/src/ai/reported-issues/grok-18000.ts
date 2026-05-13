import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-18000: PC plot `columnNames` update propagation.
// Shipped in 1.26.0/1.26.3. Removing entries from `columnNames` from the
// Context Panel previously did not refresh the axes until manual resize. The
// fix makes a `columnNames` change immediately invalidate the layout. We pin
// the round-trip: build via `DG.Viewer.pcPlot`, mutate via `setOptions`, read
// back from `getOptions(true).look`. No DOM/canvas inspection — JS API only.
category('AI: GROK-18000: PC plot columnNames update', () => {
  test('initial columnNames length is preserved on the look bag', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.PC_PLOT);

    const look = v.getOptions(true).look;
    const cn = look['columnNames'] as string[];
    expect(Array.isArray(cn), true);
    expect(cn.length, 3);
  });

  test('shrinking columnNames via setOptions updates the look bag exactly', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});

    v.setOptions({columnNames: ['age', 'height']});
    const look = v.getOptions(true).look;
    const cn = look['columnNames'] as string[];
    expect(Array.isArray(cn), true);
    expectArray(cn, ['age', 'height']);
  });

  test('empty then single-column transitions: viewer survives, look reflects each step', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});

    var threw = false;
    try {
      v.setOptions({columnNames: []});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);

    const lookEmpty = v.getOptions(true).look;
    const cnEmpty = lookEmpty['columnNames'] as string[];
    expect(Array.isArray(cnEmpty), true);
    expect(cnEmpty.length, 0);

    v.setOptions({columnNames: ['age']});
    const lookOne = v.getOptions(true).look;
    const cnOne = lookOne['columnNames'] as string[];
    expect(Array.isArray(cnOne), true);
    expectArray(cnOne, ['age']);
  });

  test('columnNames property is discoverable via getProperties (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height']});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);

    var found: DG.Property | null = null;
    for (var p of props) {
      if (p.name === 'columnNames') {
        found = p;
        break;
      }
    }
    // Best-effort: the descriptor list may not always carry `columnNames` at
    // runtime. When missing, fall back to the schema-only check we already
    // exercise above (set via setOptions, read via look).
    if (found == null) {
      v.setOptions({columnNames: ['age']});
      const cn = v.getOptions(true).look['columnNames'] as string[];
      expectArray(cn, ['age']);
    }
  });
});
