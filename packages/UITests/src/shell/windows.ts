import {category, delay, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
// import * as DG from 'datagrok-api/dg';

category('Shell: Windows', () => {
  test('ShowColumns', async () => {
    await checkSwitch('showColumns', '.d4-root.d4-column-grid');
  });
  test('ShowHelp', async () => {
    await checkSwitch('showHelp', '.grok-help-host > .grok-help');
  });
  test('ShowConsole', async () => {
    await checkSwitch('showConsole', '.d4-console-wrapper > .d4-console-header');
  });
  test('ShowProperties', async () => {
    await checkSwitch('showProperties', '.grok-prop-panel');
  });
  test('ShowRibbon', async () => {
    grok.shell.windows.simpleMode = true;
    try {
      await checkSwitch('showRibbon', '.d4-ribbon');
    } finally {
      grok.shell.windows.simpleMode = false;
    }
  });
  test('ShowSidebar', async () => {
    await checkSwitch('showSidebar', '.layout-sidebar');
  });
  test('ShowTables', async () => {
    await checkSwitch('showTables', '.grok-tables-manager');
  });
  test('ShowToolbox', async ()=> {
    await checkSwitch('showToolbox', '.d4-app-root:not(.d4-toolbox-auto)');
  });
  test('ShowVariables', async () => {
    await checkSwitch('showVariables', '.grok-view-variables.grok-variables-pane');
  });
  test('PresentationMode', async () => {
    await checkSwitch('presentationMode', '.d4-app-root.presentation');
  });
  test('SimpleMode', async () => {
    await checkSwitch('simpleMode', '.d4-app-root:not(.d4-show-menu)');
  });
}, { owner: 'aparamonov@datagrok.ai' });

function checkElementVisible(selector: string, exists: boolean = true): void {
  const e = document.body.querySelector(selector);
  if (e == undefined && exists)
    throw new Error(`Element "${selector}" not found`);
  if ((e != undefined && e instanceof HTMLElement && e.offsetParent != undefined) && !exists)
    throw new Error(`Element "${selector}" found`);
}

async function checkSwitch(sw: string, selector: string) {
  const state = (<any>grok.shell.windows)[sw];
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
    await delay(10);
  }
}
