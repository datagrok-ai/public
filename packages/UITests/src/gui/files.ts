import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, delay, test, awaitCheck} from '@datagrok-libraries/test/src/test';
import {setDialogInputValue, returnDialog, checkDialog} from './gui-utils';

const DIALOG_TITLE = 'New file share';
const CONN_NAME = 'Test Connection to S3';

category('GUI: Files', () => {
  test('demoFiles.createConnectionToS3', async () => {
    const connFilter = `dataSource in ("Dropbox","Files","Git","GitHub","GoogleCloud","S3")
and friendlyName="${CONN_NAME}"`;

    // Navigate to Files in the browse panel
    const bp = grok.shell.browsePanel;
    const filesNode = bp.mainTree.children.find((n: DG.TreeViewNode) => n.text === 'Files');
    if (filesNode == null)
      throw new Error('Files node not found in browse panel');
    filesNode.captionLabel.click();

    // Wait for the "New File Share..." button to appear
    await awaitCheck(() => {
      const btn = Array.from(document.querySelectorAll('.ui-btn-ok.ui-btn-raised'))
        .find((el) => el.textContent?.includes('New File Share'));
      return btn != null;
    }, '"New File Share..." button not found', 5000);
    const newShareBtn = Array.from(document.querySelectorAll('.ui-btn-ok.ui-btn-raised'))
      .find((el) => el.textContent?.includes('New File Share')) as HTMLElement;
    newShareBtn.click();

    await awaitCheck(() => checkDialog(DIALOG_TITLE), 'Dialog not found', 3000);
    const dialog = returnDialog(DIALOG_TITLE)!;

    // Wait for the General tab content to load asynchronously
    await awaitCheck(() => {
      try { dialog.input('Data Source'); return true; } catch (e) { return false; }
    }, 'Data Source input not loaded', 5000);
    setDialogInputValue(DIALOG_TITLE, 'Data Source', 'S3');

    // Wait for the dialog to rebuild after data source change (tabs reload asynchronously)
    await awaitCheck(() => {
      try { dialog.input('Region'); return true; } catch (e) { return false; }
    }, 'S3 connection inputs not loaded', 5000);
    setDialogInputValue(DIALOG_TITLE, 'Name', CONN_NAME);
    setDialogInputValue(DIALOG_TITLE, 'Region', 'us-east-2');
    setDialogInputValue(DIALOG_TITLE, 'Bucket', 'datagrok-data');
    setDialogInputValue(DIALOG_TITLE, 'Dir', '/packages/demo/files/demo/northwind');
    setDialogInputValue(DIALOG_TITLE, 'Access Key', '');
    setDialogInputValue(DIALOG_TITLE, 'Secret Key', '');
    await delay(300);

    const dialogRoot = dialog.root;
    const findBtn = (name: string) => Array.from(dialogRoot.querySelectorAll('.d4-command-bar .ui-btn'))
      .find((el) => el.textContent?.trim() === name) as HTMLElement | undefined;

    await awaitCheck(() => findBtn('TEST') != null, 'TEST button not found', 5000);
    findBtn('TEST')!.click();
    await awaitCheck(() => document.getElementsByClassName('d4-balloon-content').length > 0,
      'Connection test failed', 10000);
    (document.getElementsByClassName('d4-balloon-content')[0] as HTMLElement).click();

    await awaitCheck(() => findBtn('OK') != null, 'OK button not found', 3000);
    findBtn('OK')!.click();
    await delay(1000);

    const testConnection = await grok.dapi.connections.filter(connFilter).first();
    if (testConnection == undefined)
      throw new Error('Test connection was not found');
    await grok.dapi.connections.delete(testConnection);
    grok.shell.v.close();
  });
});
