import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseComparisonDefaults, entryFromDataFrame} from '../components/RunComparison/entry-extraction';

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

  test('null annotation degrades to legacy aliases', async () => {
    // JSON.parse('null') succeeds with null — must not crash on property access
    const defaults = parseComparisonDefaults({comparison: 'null', comparisonIndex: 'time'});
    expect(defaults.index, 'time');
    expect(defaults.mode == null, true);
  });

  test('missing annotations produce no defaults', async () => {
    const defaults = parseComparisonDefaults(undefined);
    expect(defaults.index == null, true);
    expect(defaults.split == null, true);
    expect(defaults.mode == null, true);
    expect(defaults.units == null, true);
  });
});

category('RunComparison: raw table entries', () => {
  const makeDf = () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'time', [1, 2, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [10, 20, 30]),
    ]);
    df.name = 'data';
    return df;
  };

  test('distinct tables with equal names and row counts get distinct ids', async () => {
    expect(entryFromDataFrame(makeDf()).id === entryFromDataFrame(makeDf()).id, false);
  });

  test('the same table maps to the same id', async () => {
    const df = makeDf();
    expect(entryFromDataFrame(df).id, entryFromDataFrame(df).id);
  });
});
