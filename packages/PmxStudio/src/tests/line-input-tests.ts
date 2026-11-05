import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {tumorModels} from '../priors';
import {LineInput, modelsFound} from '../line-input';

category('Tumor models', () => {
  test('the atlas’s models are found by class, name, tumor type or organ', async () => {
    const models = await tumorModels();
    const lines = (cls: Parameters<typeof modelsFound>[1], q: string) =>
      modelsFound(models, cls, q).flatMap((g) => g.lines);
    expect(lines(null, '').length, Object.keys(models).length, 'everything, nothing filtered');
    expect(modelsFound(models, null, '').map((g) => g.cls).join(','), 'syngeneic,pdx', 'by class, CDX empty');
    expect(lines('pdx', '').join(','), 'PDX-NSCLC,PDX-PDAC');
    expect(lines(null, 'pancrea').join(','), 'Panc02,PDX-PDAC');
    expect(lines(null, 'lung').includes('PDX-NSCLC'), true);
    expect(lines('cdx', '').length, 0);
  });

  test('the picker holds a line name, as the functions take it', async () => {
    const input = new LineInput({label: 'Tumor model'});
    expect(input.value, 'MC38');
    input.value = '4T1';
    expect(input.value, '4T1');
    input.dispose();
  });
});
