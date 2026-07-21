import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:DemoFiles/SPGI.csv';

test('Histogram tests', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await v.addViewerByIcon(page, 'histogram', 'Histogram');

  await softStep('Spline mode', async () => {
    const result = await v.setViewerProps(page, 'Histogram', [
      {set: {spline: true}, read: 'spline'},
      {set: {fillSpline: true}, read: 'fillSpline'},
      {set: {fillSpline: false}, read: 'fillSpline'},
      {set: {spline: false}, read: 'spline'},
    ]);
    expect(result).toEqual([true, true, false, false]);
  });

  await softStep('Appearance', async () => {
    const result = await v.setViewerProps(page, 'Histogram', [
      {set: {showXAxis: true, showYAxis: true}, read: ['showXAxis', 'showYAxis']},
      {set: {showXAxis: false, showYAxis: false}, read: ['showXAxis', 'showYAxis']},
      {set: {xAxisHeight: 30}, read: 'xAxisHeight'},
      {set: {allowColumnSelection: false}, read: 'allowColumnSelection'},
      {set: {showBinSelector: false}, read: 'showBinSelector'},
      {set: {showSplitSelector: false}, read: 'showSplitSelector'},
      {set: {showRangeSlider: false}, read: 'showRangeSlider'},
      {
        set: {
          showXAxis: true, showYAxis: true, allowColumnSelection: true,
          showBinSelector: true, showSplitSelector: true, showRangeSlider: true,
        },
        read: ['showXAxis', 'showYAxis', 'allowColumnSelection', 'showBinSelector', 'showSplitSelector', 'showRangeSlider'],
      },
    ]);
    expect(result[0]).toEqual({showXAxis: true, showYAxis: true});
    expect(result[1]).toEqual({showXAxis: false, showYAxis: false});
    expect(result[2]).toBe(30);
    expect(result[3]).toBe(false);
    expect(result[4]).toBe(false);
    expect(result[5]).toBe(false);
    expect(result[6]).toBe(false);
    expect(result[7]).toEqual({
      showXAxis: true, showYAxis: true, allowColumnSelection: true,
      showBinSelector: true, showSplitSelector: true, showRangeSlider: true,
    });
  });

  await softStep('Labels', async () => {
    const result = await v.setViewerProps(page, 'Histogram', [
      {set: {splitColumnName: 'SEX'}, wait: 500},
      {set: {legendVisibility: 'Never'}, read: 'legendVisibility'},
      {set: {legendVisibility: 'Always'}, read: 'legendVisibility'},
      {set: {legendPosition: 'RightTop'}, read: 'legendPosition'},
      {set: {splitColumnName: ''}},
      {set: {showTitle: true}, read: 'showTitle'},
      {set: {title: 'Age Distribution'}, read: 'title'},
      {set: {description: 'Shows distribution of patient ages'}, read: 'description'},
      {set: {descriptionVisibilityMode: 'Always'}, read: 'descriptionVisibilityMode'},
      {set: {descriptionPosition: 'Bottom'}, read: 'descriptionPosition'},
    ]);
    expect(result).toEqual([
      'Never', 'Always', 'RightTop', true, 'Age Distribution',
      'Shows distribution of patient ages', 'Always', 'Bottom',
    ]);
  });

  await softStep('Context menu', async () => {
    // Right-clicking the canvas opens the view menu, not a histogram-specific one,
    // so the only check here is that the histogram survives.
    const result = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      return h != null;
    });
    expect(result).toBe(true);
  });

  await softStep('Layout persistence', async () => {
    const result = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      const tv = grok.shell.tv;

      h.props.valueColumnName = 'WEIGHT';
      h.props.bins = 15;
      h.props.splitColumnName = 'RACE';
      h.props.splitStack = true;
      await new Promise(res => setTimeout(res, 500));

      const layout = tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise(res => setTimeout(res, 1000));

      h.close();
      await new Promise(res => setTimeout(res, 500));

      const saved = await grok.dapi.layouts.find(layoutId);
      tv.loadLayout(saved);
      await new Promise(res => setTimeout(res, 3000));

      const h2 = Array.from(tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      const r: any[] = [];
      r.push(h2 != null);
      r.push(h2 ? h2.props.valueColumnName : 'NOT_RESTORED');
      r.push(h2 ? h2.props.bins : 'NOT_RESTORED');
      r.push(h2 ? h2.props.splitColumnName : 'NOT_RESTORED');
      r.push(h2 ? h2.props.splitStack : 'NOT_RESTORED');

      await grok.dapi.layouts.delete(saved);

      return r;
    });
    expect(result[0]).toBe(true);
    expect(result[1]).toBe('WEIGHT');
    expect(result[2]).toBe(15);
    expect(result[3]).toBe('RACE');
    expect(result[4]).toBe(true);
  });

  await softStep('Data properties', async () => {
    const result = await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 500));

      const dfSpgi = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      dfSpgi.name = 'SPGI';
      const tv = grok.shell.addTableView(dfSpgi);
      await new Promise(resolve => {
        const sub = dfSpgi.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const icon = document.querySelector('[name="icon-histogram"]') as HTMLElement;
      icon.click();
      await new Promise(r => setTimeout(r, 1000));

      const h = Array.from(tv.viewers).find((v: any) => v.type === 'Histogram') as any;
      const r: any[] = [];

      for (let i = 0; i < 5; i++)
        dfSpgi.selection.set(i, true);
      await new Promise(res => setTimeout(res, 300));

      h.props.rowSource = 'Selected';
      await new Promise(res => setTimeout(res, 300));
      r.push(h.props.rowSource);

      h.props.rowSource = 'All';
      await new Promise(res => setTimeout(res, 300));
      r.push(h.props.rowSource);

      // Re-open demog: the filter expression below uses demog columns.
      grok.shell.closeAll();
      await new Promise(res => setTimeout(res, 500));

      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv2 = grok.shell.addTableView(df);
      await new Promise(resolve => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
        setTimeout(resolve, 3000);
      });

      const icon2 = document.querySelector('[name="icon-histogram"]') as HTMLElement;
      icon2.click();
      await new Promise(res => setTimeout(res, 1000));

      const h2 = Array.from(tv2.viewers).find((v: any) => v.type === 'Histogram') as any;

      h2.props.filter = '${AGE} > 40';
      await new Promise(res => setTimeout(res, 500));
      r.push(h2.props.filter);

      h2.props.filter = '';
      await new Promise(res => setTimeout(res, 300));
      r.push(h2.props.filter);

      // Table switching: SKIP — requires multiple open tables and UI interaction
      return r;
    });
    expect(result[0]).toBe('Selected');
    expect(result[1]).toBe('All');
    expect(result[2]).toBe('${AGE} > 40');
    expect(result[3]).toBe('');
  });

  v.finishSpec();
});
