import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('Grid: MultiValuesColumn', () => {
  let v: DG.TableView;
  const table = DG.DataFrame.fromCsv(`Country, Languages    
    Belgium,"Dutch
    French
    German"
    Burundi,"French
    Kirundi
    English"
    Cameroon, "English
    French"`);

    table.col('Languages')!.setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, '\n');

    before(async () => {
      v = grok.shell.addTableView(table);
    });

    test('grid.multiValuesColumn', async () => {
      const languageTags: string[] = Array.from(table.col('Languages')!.tags);

      if (languageTags[0][0] != DG.TAGS.MULTI_VALUE_SEPARATOR)
        throw 'multi-value-separator not assigned to column';
    });

    after(async () => {
      v.close();
      grok.shell.closeAll();
    });
});
