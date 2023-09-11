/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {DEFAULT_FORMATS} from '../model/const';
import {FormatConverter} from '../model/format-translation/format-converter';
import {getJsonData} from '../model/data-loading-utils/json-loader';
import {formatsToHelm} from './const';
import {_package} from '../package';

function getHelm(strand: string, format: string): string {
  return (new FormatConverter(strand, format).convertTo(DEFAULT_FORMATS.HELM));
}

category('HELM to formats', () => {
  before(async () => {
    await getJsonData();
    await _package.initMonomerLib();
  });

  for (const format of Object.keys(formatsToHelm)) {
    for (const [strand, helm] of Object.entries(formatsToHelm[format])) {
      test(`${format} to HELM`, async () => {
        const expected = helm;
        const result = getHelm(strand, format);
        expect(result, expected);
      });
    }
  }
});
