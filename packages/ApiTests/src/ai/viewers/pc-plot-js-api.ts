import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, Subscription} from 'rxjs';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: PC Plot JS API', () => {
  test('factory DG.Viewer.pcPlot returns typed DG.PcPlot with PC settings', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});
    expect(v instanceof DG.PcPlot, true);
    expect(v.type, DG.VIEWER.PC_PLOT);
    expect(v.dataFrame === df, true);
    const cols = v.props['columnNames'];
    expect(Array.isArray(cols), true);
    expect(cols.indexOf('age') >= 0, true);
    expect(cols.indexOf('height') >= 0, true);
    expect(cols.indexOf('weight') >= 0, true);
  });

  test('columnNames array round-trip via getOptions and setOptions', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height']});
    const cols = v.props['columnNames'];
    expect(Array.isArray(cols), true);
    expectArray(cols, ['age', 'height']);
    const look = v.getOptions(true).look;
    expect(Array.isArray(look['columnNames']), true);
    expectArray(look['columnNames'] as string[], ['age', 'height']);
    v.setOptions({columnNames: ['age', 'weight']});
    const cols2 = v.props['columnNames'];
    expect(Array.isArray(cols2), true);
    expectArray(cols2, ['age', 'weight']);
    const look2 = v.getOptions(true).look;
    expectArray(look2['columnNames'] as string[], ['age', 'weight']);
  });

  test('PC-specific Boolean props round-trip and getProperties surfaces columnNames', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});
    const initialNormalize = v.props['normalizeEachColumn'];
    const initialShowCurrent = v.props['showCurrentLine'];
    const initialShowAll = v.props['showAllLines'];
    v.setOptions({
      normalizeEachColumn: !initialNormalize,
      showCurrentLine: !initialShowCurrent,
      showAllLines: !initialShowAll,
    });
    expect(v.props['normalizeEachColumn'], !initialNormalize);
    expect(v.props['showCurrentLine'], !initialShowCurrent);
    expect(v.props['showAllLines'], !initialShowAll);
    const look = v.getOptions(true).look;
    expect(look['normalizeEachColumn'], !initialNormalize);
    expect(look['showCurrentLine'], !initialShowCurrent);
    expect(look['showAllLines'], !initialShowAll);
    const props = v.getProperties();
    var columnNamesProp: DG.Property | undefined;
    for (var p of props)
      if (p.name === 'columnNames') columnNamesProp = p;
    expect(columnNamesProp != null, true);
  });

  test('onLineClicked and onLineHovered are rxjs Observables', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height']}) as DG.PcPlot;
    expect(v instanceof DG.PcPlot, true);
    const subs: Subscription[] = [];
    try {
      const streams: Observable<any>[] = [v.onLineClicked, v.onLineHovered];
      for (var s of streams) {
        expect(s != null, true);
        expect(typeof s.subscribe, 'function');
        const sub = s.subscribe(() => {});
        expect(typeof sub.unsubscribe, 'function');
        subs.push(sub);
      }
    }
    finally {
      for (var sub of subs)
        sub.unsubscribe();
    }
  });

  test('view.addViewer attaches a typed DG.PcPlot', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.PC_PLOT,
        {columnNames: ['age', 'height', 'weight']}) as DG.PcPlot;
      expect(v.type, DG.VIEWER.PC_PLOT);
      expect(v instanceof DG.PcPlot, true);
      var found: DG.Viewer | undefined;
      for (var x of tv.viewers)
        if (x.type === DG.VIEWER.PC_PLOT) { found = x; break; }
      expect(found != null, true);
      expect(found instanceof DG.PcPlot, true);
    }
    finally {
      tv.close();
    }
  });
});
