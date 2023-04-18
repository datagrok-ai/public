import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {SequenceToSmilesConverter} from '../model/sequence-to-molfile-utils/sequence-to-smiles';
import {SYNTHESIZERS} from '../model/const';

function getSmiles(strand: string, format: string): string {
  return (new SequenceToSmilesConverter(strand, false, format)).convert();
}

const AXOLABS = SYNTHESIZERS.AXOLABS;

const axolabsToSmiles: {[key: string]: string} = {
  'usCfCfUfGfAf': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O',

  'usAfsusgsgsg': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',

  'UfUfUfsCfsuacg': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',

  'susususauasu': 'OP(=O)(S)OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O',

  'CfGfCfsGfsCf': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O',

  'acacacsacsac': 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O',

  'cccgggusug': 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O',

  'UfAfCfGfGfCfAfUf': 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O',

  '(invAb)cuCfuUfsc': 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O',

  '(invAb)usAfsucuCfuUfAfgcugUfgCfacususu': 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O',

  '(invAb)asacgGfuGfCfAfacucuauuca': 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O',
};

category('sequence-translator', () => {
  for (const strand of Object.keys(axolabsToSmiles)) {
    test(strand, async () => {
      const expected = axolabsToSmiles[strand];
      const result = getSmiles(strand, AXOLABS);
      expect(result, expected);
    });
  }
});
