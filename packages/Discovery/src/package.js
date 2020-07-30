/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();
let carsView;

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
   let sidebarPane1 = grok.shell.sidebar.addPane('FIRST', () => ui.div([
      ui.button('Honda', () => {
         setFilter('Honda');
      }),
      ui.button('Tesla', () => {
         setFilter('Tesla');
      })
   ]), ui.iconFA('smile'));
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

   let parser = document.createElement('a');
   parser.href = window.location;
   let pathSegments = parser.pathname.split('/');

   let v = grok.shell.addTableView(grok.data.testData("demog"));
   v.scatterPlot();
   grok.shell.v = v;

   carsView = grok.shell.addTableView(grok.data.testData("cars"));
   carsView.scatterPlot();
   carsView.basePath = '/cars'
   grok.shell.v = carsView;

   if (pathSegments.length > 4) {
      console.log(pathSegments[4]);
      setFilter(pathSegments[4]);
   }
   else
      setFilter('All');
}


function setFilter(str) {
   grok.shell.v = carsView;
   carsView.path = '/' + str;

   for (let i = 0; i < carsView.dataFrame.rowCount; i++)
      carsView.dataFrame.filter.set(i, str === 'All' || carsView.dataFrame.get('Make', i) === str);

}


//name test
//input: string s
//output: string out
//meta.progress: true
export async function countToTen(s, progress) {
   for(let i = 0; i < 10; i++) {
      if (progress != null) {
         progress.update(i * 10);
         progress.log(`Changed progress to ${i*10}`);
      }
      await new Promise(resolve => setTimeout(resolve, 1000));
   }
   return `result: ${s}`
}