/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: Discovery
//tags: app
export function startApp() {
   grok.shell.windows.showToolbox = false;
   grok.shell.windows.showHelp = false;
   grok.shell.windows.showConsole = false;

   //add top panel content
   grok.shell.topPanel.appendChild(ui.h1('Discovery'));

   //clear Datagrok native icons from sidebar
   grok.shell.sidebar.clear();

   //add custom sidebar panels
   let pane1 = grok.shell.sidebar.addPane('FIRST', () => ui.divText('A panel'), ui.iconFA('smile'));
   let pane2 = grok.shell.sidebar.addPane('SECOND', () => ui.divText('Another panel'), ui.iconFA('function'));
   grok.shell.sidebar.currentPane = pane1;

   let v = grok.shell.addTableView(grok.data.testData("demog"));
   v.scatterPlot();
}