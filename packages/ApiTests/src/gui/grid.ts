import {after, before, category, delay, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { getHTMLElementbyInnerText } from './gui-utils';
import {checkDialog} from './gui-utils';

category('GUI: Grid', () => {

  test('grid.dataSearch', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    let searchTab = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Search');
    searchTab!.click();
    
    let searchInput:HTMLInputElement | undefined;
    let input;
    for (let i=0; i<document.getElementsByClassName('ui-input-editor').length; i++) {
        input = document.getElementsByClassName('ui-input-editor')[i] as HTMLInputElement;
        if (input.placeholder == 'Search'){
            searchInput = input;
            break;
        }
    }

    await awaitCheck(() => {return searchInput != undefined});
    searchInput!.value = 'Asian'; await delay(100);
    searchInput!.dispatchEvent(new Event('input'));

    let searchOptionsBtn = document.getElementsByClassName('d4-flex-row d4-flew-nowrap d4-search')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    searchOptionsBtn.click(); 
    
    await awaitCheck(() => {return getHTMLElementbyInnerText('d4-menu-item-label', 'Filter matching') != undefined});    
    
    let filterMatchingAction = getHTMLElementbyInnerText('d4-menu-item-label', 'Filter matching'); 
    filterMatchingAction!.click();
    
    await awaitCheck(() => {return demog.filter.trueCount == 762});

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });

  test('grid.deleteRows', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);

    grok.shell.topMenu.find('Select').find('Random...').click(); 
    await awaitCheck(() => {return checkDialog('Select Random Rows')});    
  
    const input25 = Array.from(document.querySelectorAll('.d4-link-label'))
    .find((el) => el.textContent === '25%') as HTMLElement;
    input25!.click();
    await awaitCheck(() => {return demog.selection.trueCount == 250})

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await delay(200);

    let removeRowsBtn = document.getElementsByClassName('svg-remove-selected-rows')[0] as HTMLElement;
    removeRowsBtn.click(); 

    await awaitCheck(() => {return demog.rowCount == 750})

    if (demog.rowCount == 1000)
        throw 'rows are not deleted';
    
    grok.shell.closeTable(demog);
    v.close();
  });

  test('grid.deleteCols', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.o = grok.shell.t.col('disease');
    await awaitCheck(() => {return Array.from(document.querySelectorAll('.d4-accordion-title'))
    .find((el) => el.textContent == 'disease') != undefined});

    let actionsSectionOnPP = Array.from(document.getElementsByClassName('grok-entity-prop-panel')[0].querySelectorAll('.d4-accordion-pane-header'))
    .find((el) => el.textContent === 'Actions') as HTMLElement;
    if (Array.from(actionsSectionOnPP.classList).find((el) => el == 'expanded') == undefined){
      actionsSectionOnPP.click();
    }
    await awaitCheck(() => {return Array.from(actionsSectionOnPP.classList).find((el) => el == 'expanded') != undefined });
  
    let removeLinkAction = Array.from(document.querySelectorAll('.d4-link-action'))
    .find((el) => el.textContent === 'Remove') as HTMLElement;
    removeLinkAction!.click(); 
    
    await awaitCheck(() => {return demog.columns.byName('disease') == null}); 
    
    if (demog.columns.byName('disease') != null)
      throw 'disease column was not deleted'
    
    grok.shell.closeTable(demog);
    v.close();  
  }); 

  test('grid.filters', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Select').find('Random...').click(); 
    await awaitCheck(() => {return checkDialog('Select Random Rows')});   
  
    const input25 = Array.from(document.querySelectorAll('.d4-link-label'))
    .find((el) => el.textContent === '25%') as HTMLElement;
    input25!.click();
    await awaitCheck(() => {return demog.selection.trueCount == 250})

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await delay(200);

    let actionsSectionOnPP = Array.from(document.getElementsByClassName('grok-prop-panel')[0].querySelectorAll('.d4-accordion-pane-header'))
    .find((el) => el.textContent === 'Actions') as HTMLElement;
    if(Array.from(actionsSectionOnPP.classList).find((e) => e == 'expanded') == undefined)
      actionsSectionOnPP.click();  
      
    await awaitCheck(() => {return Array.from(actionsSectionOnPP.classList).find((e) => e == 'expanded') != undefined});  

    let filterLinkAction = Array.from(document.querySelectorAll('.d4-link-action'))
    .find((el) => el.textContent === 'Filter Rows') as HTMLElement;

    filterLinkAction!.click();
    await awaitCheck(() => {return demog.filter.trueCount == 250});

    if (demog.filter.trueCount != 250)
     throw 'Error in filtering'

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t)); 
  });
});
