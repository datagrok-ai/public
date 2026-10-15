import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {bulkTranslate} from '../apps/pattern/model/translator';
import {EventBus} from '../apps/pattern/model/event-bus';

category('Pattern: Bulk translate', () => {
  test('Convert twice adds a column with an unused name', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ID', ['1', '2']),
      DG.Column.fromStrings('Sense', ['', '']),
    ]);
    const eventBus = {
      getTableSelection: () => df,
      isAntisenseStrandActive: () => false,
      getSelectedStrandColumn: (strand: string) => strand === 'SS' ? 'Sense' : null,
      getSelectedIdColumn: () => 'ID',
      getPatternName: () => 'Pattern',
      getNucleotideSequences: () => ({SS: [], AS: []}),
      getTerminalModifications: () => ({SS: {'5\'': '', '3\'': ''}, AS: {'5\'': '', '3\'': ''}}),
      getPhosphorothioateLinkageFlags: () => ({SS: [false], AS: [false]}),
    } as unknown as EventBus;

    bulkTranslate(eventBus);
    bulkTranslate(eventBus);

    expect(df.columns.names().join(','), 'ID,Sense,Pattern(Sense),Pattern(Sense) (2)');
  });
});
