/* PeptiHit's campaigns, made and removed as Hit Design's are (bindings/campaign.ts), with a Helm column the app
   reads as HELM; and the HELM editor a Helm cell opens: Helm's full-screen dialog with no title, whose parts carry the
   hwe library's data-testid rather than the u2 contract, so the steps reach them by those. */
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {expect, pollMs} from '@datagrok-libraries/bdd/runtime';
import {makeCampaign} from './campaign.js';

const EDITOR = '.d4-dialog-full-screen:has([data-testid="app-root"])';

export const peptiHitCampaign = Given('the PeptiHit campaign {string} with the stages {string} and a table of:',
  async (page: Page, name: string, stages: string, table: string[][]) => {
    await makeCampaign(page, 'PeptiHit', name, stages.split(',').map((s) => s.trim()).filter((s) => s.length > 0), table);
  }, {tier: 'api', description: 'as the Hit Design one, in PeptiHit\'s folder and under its app name in the hitdesign database: the rows given (a Helm column the app tags as HELM), the template embedded with the stages given, a Chemist field, the Chemical Properties compute (MW) on the molecule the HELM converts to, the Demo File Submit preset; removed at feature end with its hitdesign rows and the V-iDs it registered, read back, and swept before it is made'});

export const helmEditorOpen = Then('the HELM editor should be open', async (page: Page) => {
  await expect(page.locator(EDITOR).filter({visible: true}), 'the HELM editor').toHaveCount(1, {timeout: pollMs(30000)});
}, {tier: 'ui', description: 'the full-screen HELM editor dialog a Helm cell opens (Helm:editMoleculeCell)'});

export const writeHelm = When('user writes the HELM {string} in the HELM editor', async (page: Page, helm: string) => {
  const editor = page.locator(EDITOR).filter({visible: true});
  await editor.locator('[data-testid="tab-helm"]').click();
  const pane = editor.locator('[data-testid="notation-pane-content"]');
  await pane.selectText();
  await page.keyboard.type(helm);
  await pane.press('Enter');
  await expect(pane, 'the HELM notation of the editor').toHaveText(helm);
  await expect(editor.locator('[data-testid="notation-pane-error"]'), 'the notation error of the editor').toHaveText('');
}, {tier: 'ui', description: 'on the editor\'s HELM tab, the notation selected and typed over, committed with Enter; done when the notation holds it and the editor reports no parse error'});

export const helmEditorOk = When('user clicks on OK in the HELM editor', async (page: Page) => {
  const editor = page.locator(EDITOR).filter({visible: true});
  await editor.locator('.d4-dialog-footer [name="button-OK"]').click();
  await expect(page.locator(EDITOR).filter({visible: true}), 'the HELM editor, closed by OK').toHaveCount(0);
}, {tier: 'ui', description: 'the editor\'s OK, which writes its HELM into the cell; done when the editor is gone'});
