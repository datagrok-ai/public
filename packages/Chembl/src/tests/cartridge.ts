import * as grok from 'datagrok-api/grok';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

category('Cartridge usage', () => {
  test('Similarity pattern', async () => {
    const df =
      await grok.data.query(`${_package.name}:patternSimilaritySearch`, {'pattern': 'c1ccc(O)cc1', 'maxRows': 1000});

    expect(df!.rowCount, 25);

    if (df != null)
      grok.shell.closeTable(df);
  });

  test('Similarity threshold pattern', async () => {
    const df = await grok.data.query(`${_package.name}:patternSimilaritySearchWithThreshold`,
      {'pattern': 'Cc1ccc2nc(N(C)CC(=O)O)sc2c1', 'threshold': 0.5});

    expect(df!.rowCount, 21);

    if (df != null)
      grok.shell.closeTable(df);
  }, {stressTest: true});

  test('Substructure pattern', async () => {
    const df =
      await grok.data.query(`${_package.name}:patternSubstructureSearch`, {'pattern': 'c1ccccc1', 'maxRows': 1000});

    expect(df!.rowCount, 1000);

    if (df != null)
      grok.shell.closeTable(df);
  }, {stressTest: true});
});
