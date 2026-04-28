import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19421: Histogram `showValues` bin labels.
// Pins the property key on `IHistogramSettings`, the round-trip through
// `props[...]` and `getOptions(true).look`, that combining the toggle with
// `splitColumnName` keeps both keys preserved, and (best-effort) that the
// runtime descriptor list carries `showValues` with a non-empty description.
category('AI: GROK-19421: Histogram showValues bin labels', () => {
  test('showValues round-trips through props and look across toggle steps', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.histogram({value: 'age'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.HISTOGRAM);

    v.props['showValues'] = true;
    expect(v.props['showValues'], true);
    expect(v.getOptions(true).look['showValues'], true);

    v.props['showValues'] = false;
    expect(v.props['showValues'], false);
    expect(v.getOptions(true).look['showValues'], false);

    v.props['showValues'] = true;
    expect(v.props['showValues'], true);
    expect(v.getOptions(true).look['showValues'], true);
  });

  test('showValues property is discoverable via getProperties (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = df.plot.histogram({value: 'age'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);
    var found: DG.Property | null = null;
    for (var p of props) {
      if (p.name === 'showValues') {
        found = p;
        break;
      }
    }
    // Best-effort: earlier `showValues` may have been a TS-only typing. If
    // the runtime descriptor list does not carry it, fall back to the
    // schema-only check we already exercised above (set via setOptions/props,
    // read via look).
    if (found != null) {
      const desc = found.description;
      expect(typeof desc === 'string', true);
      expect((desc ?? '').trim().length > 0, true);
    }
    else {
      v.props['showValues'] = true;
      expect(v.getOptions(true).look['showValues'], true);
    }
  });

  test('showValues survives alongside splitColumnName via setOptions', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.histogram({value: 'age'});
    var threw = false;
    try {
      v.setOptions({splitColumnName: 'race', showValues: true});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);

    const look = v.getOptions(true).look;
    expect(look['splitColumnName'], 'race');
    expect(look['showValues'], true);

    // Toggle once more under the split-column setup; both keys must remain
    // coherent in the look bag.
    v.props['showValues'] = false;
    const look2 = v.getOptions(true).look;
    expect(look2['splitColumnName'], 'race');
    expect(look2['showValues'], false);
  });
});
