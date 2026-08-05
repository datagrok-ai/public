import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseComparisonDefaults} from '../components/RunComparison/entry-extraction';

category('RunComparison: annotation defaults', () => {
  test('parses the comparison JSON annotation', async () => {
    const defaults = parseComparisonDefaults({
      comparison: '{"index": "time", "split": "species", "mode": "timeseries", "units": "s"}',
    });
    expect(defaults.index, 'time');
    expect(defaults.split, 'species');
    expect(defaults.mode, 'timeseries');
    expect(defaults.units, 's');
  });

  test('ignores unknown mode and units values', async () => {
    const defaults = parseComparisonDefaults({comparison: '{"mode": "spiral", "units": "weeks"}'});
    expect(defaults.mode == null, true);
    expect(defaults.units == null, true);
  });

  test('falls back to legacy aliases', async () => {
    const defaults = parseComparisonDefaults({comparisonIndex: 'time', comparisonSplit: 'species'});
    expect(defaults.index, 'time');
    expect(defaults.split, 'species');
    expect(defaults.mode == null, true);
  });

  test('comparison JSON overrides legacy aliases', async () => {
    const defaults = parseComparisonDefaults({
      comparison: '{"index": "t", "split": "kind"}',
      comparisonIndex: 'time',
      comparisonSplit: 'species',
    });
    expect(defaults.index, 't');
    expect(defaults.split, 'kind');
  });

  test('malformed JSON degrades to legacy aliases', async () => {
    const defaults = parseComparisonDefaults({comparison: '{index: time', comparisonIndex: 'time'});
    expect(defaults.index, 'time');
    expect(defaults.split == null, true);
  });

  test('missing annotations produce no defaults', async () => {
    const defaults = parseComparisonDefaults(undefined);
    expect(defaults.index == null, true);
    expect(defaults.split == null, true);
    expect(defaults.mode == null, true);
    expect(defaults.units == null, true);
  });
});
