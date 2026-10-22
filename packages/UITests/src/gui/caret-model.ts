import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

import { getHTMLElementbyInnerText } from './gui-utils';
import { after, before, category, expect, test, waitUntil, waitUntilAsync } from '@datagrok-libraries/test/src/test';

category('Modeling: Caret', () => {
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    grok.shell.addTableView(demog);
  });

  test('modeling.caret', async () => {
    grok.shell.topMenu.find('Tools').find('Predictive Modeling').find('Train Model...').click();
    await waitUntil(() => grok.shell.v.name === 'Predictive model', 'Predictive model view is not open');
    expect(grok.shell.v.name, 'Predictive model');

    function returnInputField(classList: string, inputName: string) {
      const inputFields = document.getElementsByClassName(classList);
      let input;
      for (let i = 0; i <= inputFields.length; i++) {
        input = inputFields[i].childNodes[0] as HTMLElement;
        if (input.innerText == inputName) {
          input = inputFields[i] as HTMLElement;
          break;
        }
      }
      return input;
    }

    const modelEngine = returnInputField('ui-input-choice ui-input-root',
      'Model Engine')!.childNodes[1] as HTMLInputElement;
    modelEngine.value = 'Caret';
    modelEngine.dispatchEvent(new Event('input'));

    const name = returnInputField('ui-input-text ui-input-root', 'Name')!.childNodes[1] as HTMLInputElement;
    name.value = 'Test Model Caret';
    name.dispatchEvent(new Event('input'));

    let predict = returnInputField('ui-input-root ui-input-column', 'Predict')!.childNodes[1] as HTMLElement;
    predict = predict.getElementsByClassName('d4-column-selector-column')[0] as HTMLElement;
    predict.innerText = 'sex';
    predict.dispatchEvent(new Event('input'));

    const features = returnInputField('ui-input-root ui-input-columns', 'Features')!.childNodes[1] as HTMLElement;
    features.click();
    await waitUntil(() => getHTMLElementbyInnerText('ui-btn ui-btn-ok enabled', 'OK') !== undefined,
      'Select columns dialog is not open');
    const okButtonInSelectColumn = getHTMLElementbyInnerText('ui-btn ui-btn-ok enabled', 'OK') as HTMLElement;
    okButtonInSelectColumn.click();
    await waitUntil(() => !okButtonInSelectColumn.isConnected, 'Select columns dialog is not closed');
    features.dispatchEvent(new Event('input'));

    await waitUntil(() => getHTMLElementbyInnerText('ui-btn ui-btn-ok enabled', 'Train') !== undefined,
      'Train button is not enabled');
    getHTMLElementbyInnerText('ui-btn ui-btn-ok enabled', 'Train')!.click();

    await waitUntilAsync(async () => await grok.dapi.models.filter('Test Model Caret').first() != undefined,
      'Caret Model has not been trained', 60000, 500);
  });

  after(async () => {
    await grok.dapi.models.delete(await grok.dapi.models.filter('Test Model Caret').first());
  });
}, { owner: 'oserhiienko@datagrok.ai' });
