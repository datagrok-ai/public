/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {DEFAULT_FORMATS} from '../apps/common/model/const';
import {loadJsonData} from '../apps/common/data-loader/json-loader';
import {formatsToHelm} from './const';
import {SequenceValidator} from '../apps/common/model/parsing-validation/sequence-validator';
import {getTranslatedSequences} from '../apps/translator/model/conversion-utils';
import {_package} from '../package';

function getTranslationObject(sequence: string, format: string): {[format: string]: string} {
  const indexOfInvalidChar = (new SequenceValidator(sequence)).getInvalidCodeIndex(format);
  return getTranslatedSequences(sequence, indexOfInvalidChar, format);
}

const inputs = {
  [DEFAULT_FORMATS.AXOLABS]: 'Afcgacsu',
  [DEFAULT_FORMATS.HELM]: 'RNA1{[fR](A)p.[25r](C)p.[25r](G)p.[25r](A)p.[25r](C)[sp].[25r](U)}$$$$'
};

category('Formats support', () => {
  before(async () => {
    await loadJsonData();
    await _package.initMonomerLib();
  });

  Object.entries(inputs).forEach(([format, sequence]) => {
    test(`All formats for ${format}`, async () => {
      const output = getTranslationObject(sequence, format);
      const result = Object.keys(output).length;
      // +1 due to nucleotides
      const expected = Object.keys(formatsToHelm).length + 1;
      expect(true, expected <= result);
    });
  });
});
