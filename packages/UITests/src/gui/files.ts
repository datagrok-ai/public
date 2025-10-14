import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, delay, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isDialogPresent, setDialogInputValue, isErrorBallon} from './gui-utils';

category('GUI: Files', () => {
  let testConnection: DG.DataConnection;

  test('demoFiles.createConnectionToS3', async () => {
    try {
      const openSection = document.getElementsByClassName('grok-icon fal fa-folder-open')[0] as HTMLElement;
      openSection.click();
      const filesMenuItem = Array.from(document.querySelectorAll('.d4-toggle-button'))
        .find((el) => el.textContent == 'Files') as HTMLElement;
      filesMenuItem!.click();
      await awaitCheck(() => Array.from(grok.shell.views).find((fv) => fv.type == 'files') !== undefined);
      const addBnt = document.getElementsByClassName('d4-ribbon')[0]
        .getElementsByClassName('grok-icon fal fa-plus')[0] as HTMLElement;
      addBnt.click();
      await awaitCheck(() => isDialogPresent('New file share'));
      setDialogInputValue('New file share', 'Data Source', 'S3');
      await delay(300);
      setDialogInputValue('New file share', 'Name', 'Test Connection to S3');
      setDialogInputValue('New file share', 'Region', 'us-east-2');
      setDialogInputValue('New file share', 'Bucket', 'datagrok-data');
      setDialogInputValue('New file share', 'Dir', '/packages/demo/files/demo/northwind');
      setDialogInputValue('New file share', 'Access Key', '');
      setDialogInputValue('New file share', 'Secret Key', '');
      await delay(300);
      const testBtn = Array.from(document.querySelectorAll('.ui-btn-ok'))
        .find((el) => el.textContent == 'TEST') as HTMLElement;
      testBtn.click();
      await awaitCheck((() => document.getElementsByClassName('d4-balloon-content').length > 0),
        'Connection test failed', 1500);
      isErrorBallon('"Test Connection to S3": connected successfully');
      const okBtn = Array.from(document.querySelectorAll('.ui-btn-ok'))
        .find((el) => el.textContent == 'OK') as HTMLElement;
      okBtn.click();
      await awaitCheck(() => Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent == 'Test Connection to S3') !== undefined,
      'test connection not found in tree', 1500);
      testConnection = await grok.dapi.connections.filter(`
dataSource in ("Dropbox","Files","Git","GitHub","GoogleCloud","S3")
and friendlyName="Test Connection to S3"`).first();
      if (testConnection == undefined)
        throw new Error('Test connection was not found');
    } finally {
      await grok.dapi.connections.delete(testConnection);
      testConnection = await grok.dapi.connections.filter(`
dataSource in ("Dropbox","Files","Git","GitHub","GoogleCloud","S3")
and friendlyName="Test Connection to S3"`).first();
      if (testConnection != undefined)
        throw new Error('Test connection has not been deleted');
    }
  });
});
