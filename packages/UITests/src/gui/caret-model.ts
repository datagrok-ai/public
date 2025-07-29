import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

import { getHTMLElementbyInnerText } from './gui-utils';
import { after, before, category, delay, expect, test } from '@datagrok-libraries/utils/src/test';

category('Modeling: Caret', () => {
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    grok.shell.addTableView(demog);
  });

  test('modeling.caret', async () => {
    grok.shell.topMenu.find('Tools').find('Predictive Modeling').find('Train Model...').click(); await delay(1000);
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
    modelEngine.value = 'Caret'; await delay(500);
    modelEngine.dispatchEvent(new Event('input')); await delay(500);

    const name = returnInputField('ui-input-text ui-input-root', 'Name')!.childNodes[1] as HTMLInputElement;
    name.value = 'Test Model Caret'; await delay(500);
    name.dispatchEvent(new Event('input')); await delay(500);

    let predict = returnInputField('ui-input-root ui-input-column', 'Predict')!.childNodes[1] as HTMLElement;
    predict = predict.getElementsByClassName('d4-column-selector-column')[0] as HTMLElement;
    predict.innerText = 'sex'; await delay(500);
    predict.dispatchEvent(new Event('input')); await delay(500);

    const features = returnInputField('ui-input-root ui-input-columns', 'Features')!.childNodes[1] as HTMLElement;
    features.click();
    const okButtonInSelectColumn = getHTMLElementbyInnerText('ui-btn ui-btn-ok enabled', 'OK') as HTMLElement;
    okButtonInSelectColumn.click(); await delay(500);
    features.dispatchEvent(new Event('input')); await delay(500);

    const TrainButton = getHTMLElementbyInnerText('ui-btn ui-btn-ok enabled', 'Train') as HTMLElement;
    TrainButton.click(); await delay(8000);

    if (await grok.dapi.models.filter('Test Model Caret').first() == undefined)
      throw new Error('Caret Model has not been trained');
  });

  after(async () => {
    await grok.dapi.models.delete(await grok.dapi.models.filter('Test Model Caret').first());
  });
}, { owner: 'oserhiienko@datagrok.ai' });
