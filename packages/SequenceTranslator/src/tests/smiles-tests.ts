/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {DEFAULT_FORMATS} from '../model/const';
import {getJsonData} from '../model/data-loading-utils/json-loader';
import {axolabsToSmiles} from './const';
import {_package} from '../package';
import {SequenceToMolfileConverter} from '../model/sequence-to-structure-utils/sequence-to-molfile';

function getSmiles(strand: string, format: string): string {
  const molfile = (new SequenceToMolfileConverter(strand, false, format)).convert();
  return DG.chem.convert(molfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
}

const AXOLABS = DEFAULT_FORMATS.AXOLABS;

category('Axolabs to smiles', () => {
  before(async () => {
    await getJsonData();
    await _package.initMonomerLib();
  });

  for (const strand of Object.keys(axolabsToSmiles)) {
    test(`${strand} to SMILES`, async () => {
      const expected = axolabsToSmiles[strand];
      const result = getSmiles(strand, AXOLABS);
      expect(result, expected);
    });
  }
});
