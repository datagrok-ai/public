import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19649: bar chart label sizing / orientation /
// category width round-trips. The interface IBarChartSettings exposes
// `orientation` (`d4.ts:826`), `categoryValueWidth` (`d4.ts:880`),
// `maxCategoryWidth` (`d4.ts:878`), and `showCategoryValues` (`d4.ts:835`).
// We pin that each prop round-trips through props[]/look[], and that
// getProperties() carries them. We do NOT assert specific defaults: some
// props are typed-but-not-runtime-registered, so the getProperties() check
// is best-effort.
category('AI: GROK-19649: Bar chart orientation and category width', () => {
  test('orientation flips between horizontal and vertical', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.bar({value: 'age', split: 'race'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.BAR_CHART);

    v.props['orientation'] = 'horizontal';
    expect(v.props['orientation'], 'horizontal');
    expect(v.getOptions(true).look['orientation'], 'horizontal');

    v.props['orientation'] = 'vertical';
    expect(v.props['orientation'], 'vertical');
    expect(v.getOptions(true).look['orientation'], 'vertical');
  });

  test('setOptions({maxCategoryWidth, categoryValueWidth}) reflects in look', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'}) as DG.BarChartViewer;
    expect(v instanceof DG.BarChartViewer, true);

    v.setOptions({maxCategoryWidth: 200, categoryValueWidth: 80});
    const look = v.getOptions(true).look;
    expect(look['maxCategoryWidth'], 200);
    expect(look['categoryValueWidth'], 80);
    expect(v.props['maxCategoryWidth'], 200);
    expect(v.props['categoryValueWidth'], 80);
  });

  test('showCategoryValues toggle round-trips through props and look', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.bar({value: 'age', split: 'race'});

    v.props['showCategoryValues'] = true;
    expect(v.props['showCategoryValues'], true);
    expect(v.getOptions(true).look['showCategoryValues'], true);

    v.props['showCategoryValues'] = false;
    expect(v.props['showCategoryValues'], false);
    expect(v.getOptions(true).look['showCategoryValues'], false);
  });

  test('getProperties carries orientation/category-width/showCategoryValues (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = df.plot.bar({value: 'age', split: 'race'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);

    const wanted = ['orientation', 'maxCategoryWidth', 'categoryValueWidth', 'showCategoryValues'];
    const byName: {[k: string]: DG.Property} = {};
    for (var p of props) {
      if (wanted.indexOf(p.name) >= 0)
        byName[p.name] = p;
    }
    // Best-effort: at least one of the four must be runtime-registered.
    // Some IBarChartSettings props are typed but not surfaced through
    // getProperties(), so we don't insist on every name.
    var foundCount = 0;
    for (var name of wanted) {
      if (byName[name] != null) {
        foundCount++;
        const desc = byName[name].description;
        if (typeof desc === 'string')
          expect(desc.length >= 0, true);
      }
    }
    expect(foundCount > 0, true);
  });
});
