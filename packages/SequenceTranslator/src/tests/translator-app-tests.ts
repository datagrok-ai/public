import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {after, awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';

import {DEFAULT_FORMATS} from '../apps/common/model/const';
import {DEFAULT_AXOLABS_INPUT} from '../apps/common/view/const';
import {_package} from '../package-test';
import {getHelm} from './utils';

category('Oligo Translator', () => {
  after(async () => {
    DG.Balloon.closeAll();
    grok.shell.closeAll();
  });

  test('Get SDF saves the drawn structure of a HELM input', async () => {
    const th = await _package.getTranslationHelper();
    const view: DG.ViewBase = await grok.functions.call('SequenceTranslator:oligoTranslatorApp');
    grok.shell.addView(view);
    await awaitCheck(() => $(view.root).find('textarea').length > 0, 'Sequence input did not appear', 5000);

    const textArea = $(view.root).find('textarea')[0] as HTMLTextAreaElement;
    textArea.value = getHelm(DEFAULT_AXOLABS_INPUT, DEFAULT_FORMATS.AXOLABS, th);
    textArea.dispatchEvent(new Event('input'));
    await awaitCheck(() => ($(view.root).find('select')[0] as HTMLSelectElement).value === DEFAULT_FORMATS.HELM,
      'HELM input was not detected', 5000);

    DG.Balloon.closeAll();
    $(view.root).find('button').filter((_, b) => b.textContent === 'Get SDF').trigger('click');
    expect($('.d4-balloon.warning').length, 0);
  });
});
