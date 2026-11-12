import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {DEFAULT_FORMATS} from '../apps/common/model/const';
import {getTranslatedSequences} from '../apps/translator/model/conversion-utils';
import {ITranslationHelper} from '../types';

import {_package} from '../package-test';
import {formatsToHelm} from './const';


function getTranslationObject(sequence: string, format: string, th: ITranslationHelper): { [format: string]: string } {
  const indexOfInvalidChar = th.createSequenceValidator(sequence).getInvalidCodeIndex(format);
  return getTranslatedSequences(sequence, indexOfInvalidChar, format, th);
}

const inputs = {
  [DEFAULT_FORMATS.AXOLABS]: 'Afcgacsu',
  [DEFAULT_FORMATS.HELM]: 'RNA1{[fR](A)p.[25r](C)p.[25r](G)p.[25r](A)p.[25r](C)[sp].[25r](U)}$$$$'
};

category('Formats support', () => {
  let th: ITranslationHelper;

  before(async () => {
    th = await _package.getTranslationHelper();
  });

  Object.entries(inputs).forEach(([format, sequence]) => {
    test(`All formats for ${format}`, async () => {
      const output = getTranslationObject(sequence, format, th);
      const result = Object.keys(output).length;
      // +1 due to nucleotides
      const expected = Object.keys(formatsToHelm).length + 1;
      expect(true, expected <= result);
    });
  });

  test('Translations of the Translator\'s Axolabs sample', async () => {
    const sequence = 'Afcgacsu';
    expect(th.createFormatDetector(sequence).getFormat(), DEFAULT_FORMATS.AXOLABS);
    expect(JSON.stringify(getTranslationObject(sequence, DEFAULT_FORMATS.AXOLABS, th)), JSON.stringify({
      'BioSpring': '27867*5',
      'HELM': 'RNA1{[fR](A)p.[25r](C)p.[25r](G)p.[25r](A)p.[25r](C)[sp].[25r](U)}$$$$',
      'Mermade12': 'IFGEfH',
      'Nucleotides': 'ACGACU',
    }));
  });

  test('A plain RNA HELM is translated to nucleotides, not to HELM', async () => {
    const sequence = 'RNA1{r(A)p.r(C)p.r(G)p.r(U)}$$$$';
    expect(th.createFormatDetector(sequence).getFormat(), DEFAULT_FORMATS.HELM);
    const output = getTranslationObject(sequence, DEFAULT_FORMATS.HELM, th);
    expect(output['Nucleotides'], 'ACGU');
    expect(DEFAULT_FORMATS.HELM in output, false, `${DEFAULT_FORMATS.HELM} among the translations`);
  });

  test('Codes the Translator cannot read leave no translations', async () => {
    const sequence = 'AfcgacsuQQ';
    expect(th.createFormatDetector(sequence).getFormat(), DEFAULT_FORMATS.AXOLABS);
    expect(th.createSequenceValidator(sequence).getInvalidCodeIndex(DEFAULT_FORMATS.AXOLABS), 8);
    expect(JSON.stringify(getTranslationObject(sequence, DEFAULT_FORMATS.AXOLABS, th)), '{}');
  });
});
