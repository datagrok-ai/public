import * as grok from 'datagrok-api/grok';
import {category, before, expect, test} from '@datagrok-libraries/test/src/test';
import {PackageFunctions} from '../package';

const sampleIdentifiers = [
  {
    'src_id': '3',
    'src_compound_id': 'U8Q',
  },
  {
    'src_id': '1',
    'src_compound_id': 'CHEMBL2262191',
  },
  {
    'src_id': '15',
    'src_compound_id': '656128',
  },
  {
    'src_id': '22',
    'src_compound_id': '3688392',
  },
  {
    'src_id': '28',
    'src_compound_id': 'Molport-002-053-961',
  },
];

category('identifiers', () => {
  before(async () => {
    grok.shell.closeAll();
  });


  test('identifiers', async () => {
    const res = await PackageFunctions.getCompoundsIds('HRYZQEDCGWRROX-UHFFFAOYSA-N');
    for (const expected of sampleIdentifiers) {
      const actual = res.find((r) => r.src_id === expected.src_id);
      expect(actual?.src_compound_id, expected.src_compound_id,
        `src_id ${expected.src_id}: expected compound_id ${expected.src_compound_id}, got ${actual?.src_compound_id}`);
    }
  }, {timeout: 30000});
});
