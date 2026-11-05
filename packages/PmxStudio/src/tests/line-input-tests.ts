import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {tumorModels} from '../priors';
import {LineInput, modelsFound, tumorsOf} from '../line-input';

category('Tumor models', () => {
  test('the atlas’s models are found by class, name, tumor type or organ', async () => {
    const models = await tumorModels();
    const lines = (cls: Parameters<typeof modelsFound>[1], q: string) =>
      modelsFound(models, cls, q).flatMap((g) => g.lines);
    expect(lines(null, '').length, Object.keys(models).length, 'everything, nothing filtered');
    expect(modelsFound(models, null, '').map((g) => g.cls).join(','), 'syngeneic,cdx,pdx', 'by class');
    expect(lines('cdx', '').join(','), 'MDA-MB-231,A549,SKOV3');
    expect(lines('pdx', '').join(','), 'PDX-NSCLC,PDX-PDAC');
    expect(lines(null, 'pancrea').join(','), 'Panc02,PDX-PDAC');
    expect(lines(null, 'lung').join(','), 'LLC,A549,PDX-NSCLC');
    expect(lines('cdx', 'breast').join(','), 'MDA-MB-231');
    expect(tumorsOf('MC38', models['MC38']), 'MC38 syngeneic tumors', 'the protocol names the model');
    expect(tumorsOf('PDX-PDAC', models['PDX-PDAC']).startsWith('PDX-PDAC patient-derived xenografts (PDX'), true);
    expect(tumorsOf('A549', models['A549']), 'A549 cell-line xenografts (CDX, lung)');
    expect(tumorsOf('X'), 'X tumors');
  });

  test('the picker holds a line name, as the functions take it', async () => {
    const input = new LineInput({label: 'Tumor model'});
    expect(input.value, 'MC38');
    input.value = '4T1';
    expect(input.value, '4T1');
    input.dispose();
  });
});
