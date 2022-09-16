import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, isErrorBallon} from './gui-utils';

category('Files: Demo Files Connection', () => {

  test('demoFiles.createConnectionToS3', async () => {
    let openSection = document.getElementsByClassName('grok-icon fal fa-folder-open')[0] as HTMLElement;
    openSection.click();

    let menus = document.getElementsByClassName('d4-toggle-button');
    let filesMenu:HTMLElement;
    for (let i = 0; i < menus.length; i++ ) {        
      filesMenu = menus[i] as HTMLElement;
      if (filesMenu.innerText == 'Files')
          break;
      }
    filesMenu!.click(); 

    await delay(3000);  

    let addBnt = document.getElementsByClassName('tab-handle-list-container')[0].getElementsByClassName('grok-icon fal fa-plus')[0] as HTMLElement;
    addBnt.click();      

    await delay(1000);  

    isDialogPresent('New file share');
    
    setDialogInputValue('New file share', 'Data Source', 'S3');

    await delay(1000);
    
    setDialogInputValue('New file share', 'Name', 'Test Connection to S3');
    setDialogInputValue('New file share', 'Region', 'us-east-2');
    setDialogInputValue('New file share', 'Bucket', 'datagrok-data');
    setDialogInputValue('New file share', 'Dir', '/packages/demo/files/demo/northwind');
    setDialogInputValue('New file share', 'Access Key', '');
    setDialogInputValue('New file share', 'Secret Key', '');

    await delay(500);
    
    let buttons = document.getElementsByClassName('ui-btn ui-btn-ok');
    let testBtn:HTMLElement;
    for (let i = 0; i < buttons.length; i++ ) {        
      testBtn = buttons[i] as HTMLElement;
      if (testBtn.innerText == 'TEST')
          break;
      }
    testBtn!.click(); await delay(1000);

    isErrorBallon('"Test Connection to S3": connected successfully');
    
    buttons = document.getElementsByClassName('ui-btn ui-btn-ok');
    let okBtn:HTMLElement;
    for (let i = 0; i < buttons.length; i++ ) {        
      okBtn = buttons[i] as HTMLElement;
      if (okBtn.innerText == 'OK')
          break;
      }
    okBtn!.click(); await delay(5000);
   
    let testConnection = await grok.dapi.connections.filter('dataSource in ("Dropbox","Files","Git","GitHub","GoogleCloud","S3") and friendlyName="Test Connection to S3"').first();

    if (testConnection == undefined)
      throw 'Test connection was not found'
  
    await grok.dapi.connections.delete(testConnection);

    testConnection = await grok.dapi.connections.filter('dataSource in ("Dropbox","Files","Git","GitHub","GoogleCloud","S3") and friendlyName="Test Connection to S3"').first();
    if (testConnection == undefined)
      throw 'Test connection has not been deleted'  
    }); 
});
