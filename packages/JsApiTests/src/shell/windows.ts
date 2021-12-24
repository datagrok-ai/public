import {after, before, category, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let v: DG.TableView;

category('Shell', () => {
  before(async () => {
    v = grok.shell.addTableView(grok.data.demo.demog());
  });

  test('Windows - ShowColumns', async () => {
    await checkSwitch('showColumns', '.d4-root.d4-column-grid');
  })
  test('Windows - ShowHelp', async () => {
    await checkSwitch('showHelp', '.grok-help-host > .grok-help');
  });
  test('Windows - ShowConsole', async () => {
    await checkSwitch('showConsole', '.d4-console-wrapper > .d4-console-header');
  });
  test('Windows - ShowProperties', async () => {
    await checkSwitch('showProperties', '.grok-prop-panel');
  });
  test('Windows - ShowRibbon', async () => {
    await checkSwitch('showRibbon', '.d4-ribbon');
  });
  test('Windows - ShowSidebar', async () => {
    await checkSwitch('showSidebar', '.layout-sidebar');
  });
  test('Windows - ShowTables', async () => {
    await checkSwitch('showTables', '.grok-tables-manager');
  });
  test('Windows - ShowToolbox', async ()=> {
    await checkSwitch('showToolbox', '.d4-app-root:not(.d4-toolbox-hidden)');
  });
  test('Windows - ShowVariables', async () => {
    await checkSwitch('showVariables', '.grok-view-variables.grok-variables-pane');
  });
  test('Windows - PresentationMode', async () => {
    await checkSwitch('presentationMode', '.d4-app-root.presentation');
  });
  test('Windows - SimpleMode', async () => {
    await checkSwitch('simpleMode', '.d4-app-root:not(.d4-show-menu)');
  });

  after(async () => {
    v.close();
    grok.shell.closeTable(v.table!);
  });

});

function checkElementVisible(selector: string, exists: boolean = true): void {
  let e = document.body.querySelector(selector);
  if (e == undefined && exists)
    throw `Element "${selector}" not found`;
  if ((e != undefined && e instanceof HTMLElement && e.offsetParent != undefined) && !exists)
    throw `Element "${selector}" found`;
}

async function checkSwitch(sw: string, selector: string) {

  let state = (<any>grok.shell.windows)[sw];
  try {
    (<any>grok.shell.windows)[sw] = true;
    checkElementVisible(selector);
    (<any>grok.shell.windows)[sw] = false;
    checkElementVisible(selector, false);
    (<any>grok.shell.windows)[sw] = true;
    checkElementVisible(selector);
  } catch (x) {
    throw x;
  } finally {
    (<any>grok.shell.windows)[sw] = state;
    await new Promise(r => setTimeout(r, 10));
  }
}
