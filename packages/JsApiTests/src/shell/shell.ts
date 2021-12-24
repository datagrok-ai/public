import {after, before, category, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let v: DG.View;

category('Shell', () => {
  before(async () => {
    v = grok.shell.addTableView(grok.data.demo.demog());
    grok.shell.windows.showColumns;
    grok.shell.windows.showHelp;
    grok.shell.windows.showConsole;
    grok.shell.windows.showProperties;
    grok.shell.windows.showRibbon;
    grok.shell.windows.showSidebar;
    grok.shell.windows.showTables;
    grok.shell.windows.showToolbox;
    grok.shell.windows.showVariables;
    grok.shell.windows.presentationMode;
    grok.shell.windows.simpleMode;
  });


  test('Windows - ShowColumns', async () => {
    checkSwitch((value) => grok.shell.windows.showColumns = value, '.d4-root.d4-column-grid');
  })
  test('Windows - ShowHelp', async ()=> {
    checkSwitch((value) => grok.shell.windows.showHelp = value, '.grok-help-host > .grok-help');
  });
  test('Windows - ShowConsole', async ()=> {
    checkSwitch((value) => grok.shell.windows.showConsole = value, '.d4-console-wrapper .d4-console-header');
  });
  test('Windows - ShowProperties', async ()=> {
    checkSwitch((value) => grok.shell.windows.showProperties = value, '.grok-prop-panel');
  });
  test('Windows - ShowRibbon', async ()=> {
    checkSwitch((value) => grok.shell.windows.showRibbon = value, '.d4-ribbon');
  });

  after(async () => {
    v.close();
  });

});

function checkElementVisible(selector: string, exists: boolean = true):void {
  let e = document.body.querySelector(selector);
  if (e == undefined && exists)
    throw `Element "${selector}" not found`;
  if ((e != undefined && e instanceof HTMLElement && e.offsetParent != undefined) && !exists)
    throw `Element "${selector}" found`;
}

function checkSwitch(sw: (arg0: boolean) => void, selector: string) {
  sw(true);
  checkElementVisible(selector);
  sw(false);
  checkElementVisible(selector, false);
}
