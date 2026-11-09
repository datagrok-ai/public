import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {MAX_SLOT_CATEGORIES, resolveSuggestion} from '../ai/prompt-suggestions';

const GROUPING_SUGGESTIONS = [
  {
    label: 'Scatter plot of {x.isNumerical} vs {y.isNumerical}, colored by {c.isCategorical}',
    prompt: 'Add a scatter plot of {x.isNumerical} vs {y.isNumerical}, colored by {c.isCategorical}',
    prefer: {c: 'name~="stage|class|status|group|type|series"'},
  },
  {prompt: 'Add a box plot of {n.isNumerical} grouped by {c.isCategorical}'},
  {prompt: 'Add a bar chart of row counts per {c.isCategorical}'},
  {prompt: 'Add a trellis plot of {n.isNumerical} split by {c.isCategorical}'},
  {
    label: 'Pivot {n.isNumerical} by {c.isCategorical}',
    prompt: 'Pivot the table — group by {c.isCategorical} and aggregate {n.isNumerical}',
  },
];

category('Prompt suggestions', () => {
  test('demog: categorical slots skip USUBJID', async () => {
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    for (const s of GROUPING_SUGGESTIONS) {
      const r = resolveSuggestion(s, df);
      expect(r !== null, true, s.prompt);
      expect(r!.prompt!.includes('USUBJID'), false, r!.prompt);
      expect(r!.label?.includes('USUBJID') ?? false, false, r!.label);
    }
  });

  test('ID-only string column: grouping suggestions dropped', async () => {
    const n = 100;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', Array.from({length: n}, (_, i) => `S-${i}`)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', Array.from({length: n}, (_, i) => i)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', Array.from({length: n}, (_, i) => i * 2)),
    ]);
    for (const s of GROUPING_SUGGESTIONS)
      expect(resolveSuggestion(s, df), null, s.prompt);
    expect(resolveSuggestion({prompt: 'Show a histogram of {n.isNumerical}'}, df) !== null, true);
  });

  test('Too many or a single category: column skipped', async () => {
    const n = 1000;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('many', Array.from({length: n}, (_, i) => `c${i % (MAX_SLOT_CATEGORIES + 1)}`)),
      DG.Column.fromStrings('constant', Array.from({length: n}, () => 'x')),
      DG.Column.fromStrings('arm', Array.from({length: n}, (_, i) => i % 2 ? 'A' : 'B')),
    ]);
    expect(resolveSuggestion({prompt: 'Count per {c.isCategorical}'}, df)?.prompt, 'Count per arm');
  });
});
