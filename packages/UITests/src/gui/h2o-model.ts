import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, getHTMLElementbyInnerText} from './gui-utils';

category('Modeling: H20', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('modeling.h2o', async () => {
    grok.shell.topMenu.find('Tools').find('Predictive Modeling').find('Train Model...').click(); await delay(1000);
    expect(grok.shell.v.name, 'Predictive model');
    
    function returnInputField(classList:string, inputName:string){
        let inputFields = document.getElementsByClassName(classList) 
        let input;
        for(let i = 0; i <= inputFields.length; i++){
            input = inputFields[i].childNodes[0] as HTMLElement;
            if (input.innerText == inputName){
                input = inputFields[i] as HTMLElement;
                break;
            }
        } 
        return input;
    }

    let modelEngine = returnInputField('ui-input-choice ui-input-root','Model Engine')!.childNodes[1] as HTMLInputElement;
    modelEngine.value = 'H2O'; await delay(500);
    modelEngine.dispatchEvent(new Event('input')); await delay(500);

    let name = returnInputField('ui-input-text ui-input-root','Name')!.childNodes[1] as HTMLInputElement;
    name.value = 'Test Model H2O'; await delay(500);
    name.dispatchEvent(new Event('input')); await delay(500);

    let predict = returnInputField('ui-input-root ui-input-column','Predict')!.childNodes[1] as HTMLElement;
    predict = predict.getElementsByClassName('d4-column-selector-column')[0] as HTMLElement;
    predict.innerText = 'sex'; await delay(500);
    predict.dispatchEvent(new Event('input')); await delay(500);

    let features = returnInputField('ui-input-root ui-input-columns','Features')!.childNodes[1] as HTMLElement;
    features.click();
    let okButtonInSelectColumn = getHTMLElementbyInnerText('ui-btn ui-btn-ok enabled', 'OK') as HTMLElement;
    okButtonInSelectColumn.click(); await delay(500);
    features.dispatchEvent(new Event('input')); await delay(500);
    

  });

  after(async () => {
    grok.shell.closeAll();
  });
});
