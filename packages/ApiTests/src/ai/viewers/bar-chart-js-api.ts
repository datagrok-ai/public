import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {take} from 'rxjs/operators';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: BarChart JS API', () => {
  test('factory typed via DG.Viewer.barChart', async () => {
    const df = grok.data.demo.demog(100);
    const v = DG.Viewer.barChart(df, {value: 'age', valueAggrType: 'avg', split: 'race'});
    expect(v.type, DG.VIEWER.BAR_CHART);
    expect(v instanceof DG.BarChartViewer, true);
    expect(v.props['valueColumnName'], 'age');
    expect(v.props['valueAggrType'], 'avg');
    expect(v.props['splitColumnName'], 'race');
    const look = v.getOptions(true).look;
    expect(look['valueColumnName'], 'age');
    expect(look['valueAggrType'], 'avg');
    expect(look['splitColumnName'], 'race');
  });

  test('factory via df.plot.bar shorthand', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.bar({value: 'height', split: 'race', stack: 'sex'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.BAR_CHART);
    expect(v.dataFrame === df, true);
    expect(v.props['valueColumnName'], 'height');
    expect(v.props['splitColumnName'], 'race');
    expect(v.props['stackColumnName'], 'sex');
  });

  test('resetView does not throw on attached BarChart', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.BAR_CHART, {value: 'age', split: 'race'}) as DG.BarChartViewer;
      expect(v.type, DG.VIEWER.BAR_CHART);
      expect(v instanceof DG.BarChartViewer, true);
      var threw = false;
      var result: any = 'sentinel';
      try {
        result = v.resetView();
      }
      catch (_e) {
        threw = true;
      }
      expect(threw, false);
      expect(result === undefined, true);
    }
    finally {
      tv.close();
    }
  });

  test('onCategoryClicked and onCategoryHovered are rxjs Observables', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'}) as DG.BarChartViewer;
    expect(v instanceof DG.BarChartViewer, true);
    const subs: Subscription[] = [];
    try {
      const streams = [v.onCategoryClicked, v.onCategoryHovered];
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

  test('onResetView round-trip via resetView', async () => {
    const df = grok.data.demo.demog(50);
    const tv = grok.shell.addTableView(df);
    let sub: Subscription | undefined;
    try {
      const v = tv.addViewer(DG.VIEWER.BAR_CHART, {value: 'age', split: 'race'}) as DG.BarChartViewer;
      expect(v instanceof DG.BarChartViewer, true);
      const fired: Promise<boolean> = new Promise((resolve) => {
        sub = v.onResetView.pipe(take(1)).subscribe(() => resolve(true));
      });
      const timeout: Promise<boolean> = (async () => {
        await DG.delay(1500);
        return false;
      })();
      v.resetView();
      const ok = await Promise.race([fired, timeout]);
      expect(ok, true);
    }
    finally {
      if (sub != null)
        sub.unsubscribe();
      tv.close();
    }
  });
});
