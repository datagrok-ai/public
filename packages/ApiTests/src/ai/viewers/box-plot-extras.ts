import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, settled, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectRoundTrip, withTableView} from '../helpers';

// BoxPlot look-prop round-trips: plotStyle, violin bins, per-statistic togglers, binColorAggrType, axisType, opacity.
category('AI: Viewers: BoxPlot extras', () => {
  const v = (): DG.BoxPlot => DG.Viewer.boxPlot(demog(), {value: 'age', category1: 'race'});

  test('plotStyle box/violin round-trip + getProperties choices', async () => {
    const c = v();
    expectRoundTrip(c, {plotStyle: 'violin'});
    expectRoundTrip(c, {plotStyle: 'box'});
    expectChoices(c, 'plotStyle', ['box', 'violin']);
  });

  test('setMarkerOpacity sets the own markerOpacity', async () => {
    const c = v();
    c.setMarkerOpacity(60);
    expect(c.props.markerOpacity, 60);
  });

  test('violin bins int round-trip within 20-1000', async () => {
    const c = v();
    c.setOptions({plotStyle: 'violin'});
    for (const n of [100, 20, 1000])
      expectRoundTrip(c, {bins: n});
  });

  test('per-statistic toggler bools combined round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {
      showStatistics: true,
      showMin: false, showMax: false, showAvg: false,
      showMed: false, showStdev: true, showQ1: true, showQ3: true,
    });
    expectRoundTrip(c, {
      showMin: true, showMax: true, showAvg: true,
      showMed: true, showStdev: false, showQ1: false, showQ3: false,
    });
  });

  test('binColorAggrType round-trip + getProperties choices', async () => {
    const c = v();
    expectRoundTrip(c, {binColorColumnName: 'height', binColorAggrType: 'med'});
    expectChoices(c, 'binColorAggrType', ['avg', 'med']);
  });

  test('axisType linear/logarithmic round-trip + getProperties choices', async () => {
    const c = v();
    expectRoundTrip(c, {axisType: 'logarithmic'});
    expectRoundTrip(c, {axisType: 'linear'});
    expectChoices(c, 'axisType', ['linear', 'logarithmic']);
  });

  test('markerOpacity 0..100 boundary round-trip', async () => {
    const c = v();
    for (const n of [0, 50, 100])
      expectRoundTrip(c, {markerOpacity: n});
  });

  test('showInsideValues/showOutsideValues/showPValue/showEmptyCategories combined bool round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {
      showInsideValues: false, showOutsideValues: false, showPValue: false, showEmptyCategories: false,
    });
    expectRoundTrip(c, {
      showInsideValues: true, showOutsideValues: true, showPValue: true, showEmptyCategories: true,
    });
  });

  test('violin over a double column holding NaN draws and takes a hover without an error (core#960)', async () => {
    const t = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, NaN, 3, 4, 5, 6, 7, 8]),
      DG.Column.fromStrings('c', ['a', 'a', 'a', 'a', 'b', 'b', 'b', 'b']),
    ]);
    expect(t.col('x')!.isNone(1), false);
    grok.shell.clearLastError();
    await withTableView(t, async (tv) => {
      const v = DG.Viewer.boxPlot(t, {valueColumnName: 'x', category1ColumnName: 'c', plotStyle: 'violin'});
      tv.addViewer(v);
      await settled(v);
      const overlay = v.root.querySelector('canvas[name="overlay"]')!;
      const r = overlay.getBoundingClientRect();
      overlay.dispatchEvent(new MouseEvent('mousemove',
        {bubbles: true, clientX: r.left + r.width / 4, clientY: r.top + r.height / 2}));
      await settled(v);
      expect(await grok.shell.lastError ?? '', '');
    });
  });

  test('a pointer over a box plot shrunk under Auto Layout logs no error (core#960)', async () => {
    grok.shell.clearLastError();
    const v = DG.Viewer.boxPlot(demog(100), {valueColumnName: 'age', category1ColumnName: 'sex',
      autoLayout: true, showStatistics: true});
    const host = ui.div([v.root], {style: {width: '170px', height: '150px'}});
    document.body.append(host);
    try {
      await settled(v);
      expect(v.getWidgetStatus().hitAreas?.['stats'] == null, true, 'the shrunk plot still lays out its statistics');
      const overlay = v.root.querySelector('canvas[name="overlay"]')!;
      const r = overlay.getBoundingClientRect();
      for (const share of [0.25, 0.5, 0.75]) {
        overlay.dispatchEvent(new MouseEvent('mousemove',
          {bubbles: true, clientX: r.left + r.width * share, clientY: r.top + r.height / 2}));
        await settled(v);
      }
      expect(await grok.shell.lastError ?? '', '');
    }
    finally {
      v.detach();
      host.remove();
    }
  });
}, {owner: 'agolovko@datagrok.ai'});
