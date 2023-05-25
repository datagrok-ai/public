import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


category('Settings', () => {
  let gss: DG.Settings & any;
  let views: DG.View[];

  before(async () => {
    gss = grok.shell.settings;
    views = [];
  });

  test('Dev Settings', async () => {
    expect(gss.loadDefaultsOnStart, true);
    expect(gss.webpackDevUrl, null);
    expect(gss.cvmUrl.startsWith('https://cvm'), true);
    expect(gss.cvmSplit, true);
    expect(gss.apiUrl.endsWith('/api'), true);
    expect(gss.helpBaseUrl, '');
    expect(gss.jupyterNotebook, `${gss.cvmUrl}/notebook`);
    expect(typeof gss.jupyterGatewayToken, 'string');
    expect(typeof gss.jupyterNotebookToken, 'string');
  });

  test('Beta Settings', async () => {
    expect(gss.enableBetaViewers, false);
    expect(gss.saveProjectWithViewLayout, false);
    try {
      gss.saveProjectWithViewLayout = true;
      expect(gss.saveProjectWithViewLayout, true);
    } finally {
      gss.saveProjectWithViewLayout = false;
    }
  });

  test('Windows Settings', async () => {
    expect(gss.showMenu, !grok.shell.windows.simpleMode,
      `showMenu ${gss.showMenu} != ${!grok.shell.windows.simpleMode}`);
    expect(gss.showTables, grok.shell.windows.showTables,
      `showTables ${gss.showTables} != ${!grok.shell.windows.showTables}`);
    expect(gss.showColumns, grok.shell.windows.showColumns,
      `showColumns ${gss.showColumns} != ${!grok.shell.windows.showColumns}`);
    expect(gss.showProperties, grok.shell.windows.showProperties,
      `showProperties ${gss.showProperties} != ${!grok.shell.windows.showProperties}`);
    expect(gss.showToolbox, grok.shell.windows.showToolbox,
      `showToolbox ${gss.showToolbox} != ${!grok.shell.windows.showToolbox}`);
    expect(gss.showStatusBar, grok.shell.windows.showStatusBar,
      `showStatusBar ${gss.showStatusBar} != ${!grok.shell.windows.showStatusBar}`);
    expect(gss.showVariables, grok.shell.windows.showVariables,
      `showVariables ${gss.showVariables} != ${!grok.shell.windows.showVariables}`);
    expect(gss.showConsole, grok.shell.windows.showConsole,
      `showConsole ${gss.showConsole} != ${!grok.shell.windows.showConsole}`);
    expect(gss.showHelp, grok.shell.windows.showHelp,
      `showHelp ${gss.showHelp} != ${!grok.shell.windows.showHelp}`);
  });

  after(async () => {
    views.forEach((v) => v.close());
    grok.shell.closeAll();
  });
});
