import {after, before, category, delay, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isViewerPresent, waitForElement} from '../gui-utils';


category('Viewers: Form', () => {
  let tv: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    tv = grok.shell.addTableView(demog);
  });

  
  test('Form.visual', async () => {
    const formIcon = document.getElementsByClassName('svg-form')[0] as HTMLElement;
    await formIcon.click();
    isViewerPresent(Array.from(tv.viewers), 'Form');

    const form = document.querySelector('[name=viewer-Form]') as HTMLElement;
    const left = form.getElementsByClassName('grok-icon fal fa-chevron-left')[0] as HTMLElement; 
    const right = form.getElementsByClassName('grok-icon fal fa-chevron-right')[0] as HTMLElement;
    const select = form.getElementsByClassName('grok-icon fal fa-square')[0] as HTMLElement;
    const edit = form.getElementsByClassName('grok-icon fal fa-edit')[0] as HTMLElement;
    const design = form.getElementsByClassName('grok-icon fal fa-object-ungroup')[0] as HTMLElement; 
    for (let i = 0; i < 10; i++) {
      right.click();
      if (i % 2 === 0) select.click();
    }
    for (let i = 0; i < 13; i++) left.click();
    
    const selected = tv.dataFrame.selection.getSelectedIndexes().length;
    expect(selected, 5);
    
    edit.click();
    Array.from(form.getElementsByClassName('d4-host-element-panel')).forEach((el, i) => {
      if (i % 2 !== 0 && i < 10)
      //@ts-ignore
      el.firstChild.value = 1;
        const  event = new Event('input', {
          bubbles: true,
          cancelable: true,
        });
        el.firstChild?.dispatchEvent(event);
    });
    
    const arr = Array.from(tv.dataFrame.currentRow.cells).map(c => c.value).slice(0, 5);
    expectArray(arr, [1, '1', '1', 1, '1']);
    
    design.click();
    
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await waitForElement('.property-grid-item-editor-checkbox', 'cannot find property panel');
      
    const checkboxes = Array.from(document.querySelectorAll('.property-grid-item-editor-checkbox'));
    for (let cb of checkboxes) {
      (cb as HTMLElement).click();
      await delay(10);
    }
  });


  after(async () => {
    tv.close();
    grok.shell.closeTable(demog);
  }); 
});
