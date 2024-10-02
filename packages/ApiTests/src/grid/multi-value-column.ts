import {category, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Grid: MultiValuesColumn', () => {
  const table = DG.DataFrame.fromCsv(`Country, Languages    
    Belgium,"Dutch
    French
    German"
    Burundi,"French
    Kirundi
    English"
    Cameroon, "English
    French"`);
  table.col('Languages')!.meta.multiValueSeparator = '\n';

  test('grid.multiValuesColumn', async () => {
    grok.shell.addTableView(table);
    const languageTags: string[] = Array.from(table.col('Languages')!.tags);

    if (languageTags[0][0] != DG.TAGS.MULTI_VALUE_SEPARATOR)
      throw new Error('multi-value-separator not assigned to column');
  });
});
