import * as grok from 'datagrok-api/grok';
import {category, before, after, expect, test} from '@datagrok-libraries/utils/src/test';
import {getCompoundsIds} from '../package';

const sampleIdentifiers = [
  {
    'src_id': '1',
    'src_compound_id': 'CHEMBL2262191',
  },
  {
    'src_id': '10',
    'src_compound_id': '12602277',
  },
  {
    'src_id': '15',
    'src_compound_id': 'SCHEMBL656128',
  },
  {
    'src_id': '3',
    'src_compound_id': 'U8Q',
  },
  {
    'src_id': '29',
    'src_compound_id': 'J1.445.979K',
  },
  {
    'src_id': '9',
    'src_compound_id': 'ZINC000002579331',
  },
  {
    'src_id': '22',
    'src_compound_id': '3688392',
  },
];

category('identifiers', () => {
  before(async () => {
    grok.shell.closeAll();
  });


  test('identifiers', async () => {
    const res: any[] | null = await getCompoundsIds('HRYZQEDCGWRROX-UHFFFAOYSA-N');
    expect(res?.length, 7, `Expected 7 identifiers, got ${res?.length}`);
    for (let i = 0; i < res!.length; i++) {
      expect(res![i].src_id, sampleIdentifiers[i].src_id);
      expect(res![i].src_compound_id, sampleIdentifiers[i].src_compound_id);
    }
  }, {timeout: 30000});


  after(async () => {
    grok.shell.closeAll();
  });
});
