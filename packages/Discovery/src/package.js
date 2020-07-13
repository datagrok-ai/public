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
   let sidebarPane1 = grok.shell.sidebar.addPane('FIRST', () => ui.divText('A panel'), ui.iconFA('smile'));
   let sidebarPane2 = grok.shell.sidebar.addPane('SECOND', () => ui.divText('Another panel'), ui.iconFA('function'));
   grok.shell.sidebar.currentPane = sidebarPane1;

   let tabs = document.createElement('div');
   tabs.classList.add('discovery-tabs');
   grok.shell.bottomPanel.appendChild(tabs);

   let pages = {};

   grok.events.onViewAdded.subscribe((view) => {
      // add page to tab control
      let page = document.createElement('div');
      page.innerText = view.name;
      page.onclick = () => grok.shell.v = view;
      tabs.appendChild(page);
      pages[view.name] = page;
   });

   grok.events.onViewRemoved.subscribe((view) => {
      pages[grok.shell.v.name].remove();
   });

   grok.events.onCurrentViewChanged.subscribe((_) => {
      [...tabs.children].forEach((p) => {
         if (p === pages[grok.shell.v.name])
            p.classList.add('current');
         else
           p.classList.remove('current');
      });
   });

   let v = grok.shell.addTableView(grok.data.testData("demog"));
   v.scatterPlot();
   grok.shell.v = v;

   let v2 = grok.shell.addTableView(grok.data.testData("cars"));
   v2.scatterPlot();
   grok.shell.v = v2;
}
