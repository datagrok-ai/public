import {after, before, category, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Shell', () => {
  before(async () => {
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
    checkSwitch((value) => grok.shell.windows.showColumns = value, '.layout-dockarea > .view-tabs > .splitter-container-row > .panel-base.splitter-container-horizontal > .panel-content > .d4-root.d4-column-grid');
  })
  test('Windows - ShowHelp', async ()=> {
    checkSwitch((value) => grok.shell.windows.showHelp = value, '.layout-dockarea > .view-tabs > .splitter-container-row > .panel-base.splitter-container-horizontal > .panel-content > .d4-root.d4-column-grid');
  });

  after(async () => {
    console.log('hi');
  });

});

function checkElementVisible(selector: string, exists: boolean = true):void {
  let e = document.body.querySelector(selector);
  if (e == undefined && exists)
    throw `Element "${selector}" not found`;
  if (e != undefined && !exists)
    throw `Element "${selector}" found`;
}

function checkSwitch(sw: (arg0: boolean) => void, selector: string) {
  sw(true);
  checkElementVisible(selector);
  sw(false);
  checkElementVisible(selector, false);
}
